"""Part 3 - Assembly pipeline: reads -> k-mers -> De Bruijn graph -> clean ->
Eulerian path (or unitig fallback) -> contig(s).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from .euler import (EulerCheck, check_euler, enumerate_euler_sequences, hierholzer,
                    path_to_sequence, unitigs)
from .graph import DeBruijnGraph
from .kmers import count_kmers, total_kmers
from .sequence import validate_dna

AMBIGUITY_CHECK_MAX_EDGES = 400


@dataclass
class AssemblyResult:
    k: int
    reads: list
    freq: dict
    graph: DeBruijnGraph
    euler: EulerCheck
    path: list = None
    contigs: list = field(default_factory=list)
    method: str = "empty"          # eulerian_path | eulerian_cycle | unitigs | empty
    ambiguous: bool = None         # True / False / None (not checked)
    alternatives: list = field(default_factory=list)
    cleaning: dict = field(default_factory=dict)
    timings: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    @property
    def success(self) -> bool:
        return self.method.startswith("eulerian")

    @property
    def sequence(self):
        """The single assembled sequence, or None if assembly fell back to unitigs."""
        return self.contigs[0] if self.success else None

    @property
    def longest_contig(self) -> str:
        return max(self.contigs, key=len) if self.contigs else ""


def _as_strings(reads) -> list:
    out = []
    for i, r in enumerate(reads, start=1):
        seq = getattr(r, "sequence", r)
        out.append(validate_dna(seq, getattr(r, "name", f"read{i}")))
    return out


def assemble(reads, k: int, min_coverage: int = 1, clean: bool = False,
             max_tip_length: int = None, check_ambiguity: bool = True,
             multigraph: bool = False) -> AssemblyResult:
    """Run the whole assembly half of the project.

    reads          list of DNA strings or SequenceRecords
    k              k-mer length (>= 2)
    min_coverage   drop k-mers seen fewer times (sequencing-error filter)
    clean          additionally remove tips and pop bubbles
    multigraph     one edge per k-mer occurrence (only for a single known sequence)
    """
    if k < 2:
        raise ValueError("k must be >= 2 (vertices are (k-1)-mers)")
    if min_coverage < 1:
        raise ValueError("min_coverage must be >= 1")
    timings, warnings = {}, []

    t0 = time.perf_counter()
    reads = _as_strings(reads)
    if not reads:
        raise ValueError("no reads given")
    short = [r for r in reads if len(r) < k]
    if short:
        warnings.append(f"{len(short)} read(s) shorter than k={k} contribute no k-mers")
    freq = count_kmers(reads, k)
    timings["kmer_counting"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    graph = DeBruijnGraph.from_kmer_counts(freq, k, min_coverage, multigraph)
    timings["graph_construction"] = time.perf_counter() - t0

    cleaning = {
        "total_kmers": total_kmers(reads, k),
        "distinct_kmers": len(freq),
        "low_coverage_removed": sum(1 for c in freq.values() if c < min_coverage),
        "tip_edges_removed": 0,
        "bubble_edges_removed": 0,
    }
    if clean:
        t0 = time.perf_counter()
        cleaning["tip_edges_removed"] = graph.remove_tips(max_tip_length)
        cleaning["bubble_edges_removed"] = graph.pop_bubbles(max_tip_length)
        cleaning["tip_edges_removed"] += graph.remove_tips(max_tip_length)
        timings["graph_cleaning"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    euler = check_euler(graph)
    path = hierholzer(graph) if euler.ok else None
    timings["eulerian_traversal"] = time.perf_counter() - t0

    result = AssemblyResult(k, reads, freq, graph, euler, path,
                            cleaning=cleaning, timings=timings, warnings=warnings)
    if path:
        result.contigs = [path_to_sequence(path)]
        result.method = "eulerian_cycle" if euler.kind == "cycle" else "eulerian_path"
        if check_ambiguity:
            if graph.num_edges() <= AMBIGUITY_CHECK_MAX_EDGES:
                t0 = time.perf_counter()
                seqs, exhaustive = enumerate_euler_sequences(graph, limit=5, max_steps=50_000)
                timings["ambiguity_check"] = time.perf_counter() - t0
                if len(seqs) > 1:
                    result.ambiguous = True
                    result.alternatives = sorted(s for s in seqs if s != result.contigs[0])
                    warnings.append(
                        f"{'at least ' if not exhaustive else ''}{len(seqs)} different sequences "
                        "use exactly these k-mers (repeat ambiguity); try a larger k")
                elif exhaustive:
                    result.ambiguous = False
            else:
                warnings.append("ambiguity check skipped (graph too large for brute force)")
        if euler.kind == "cycle":
            warnings.append("graph is an Eulerian cycle: first and last (k-1)-mers coincide, "
                            "so the start point (rotation) is a convention")
    elif graph.num_edges() > 0:
        t0 = time.perf_counter()
        result.contigs = unitigs(graph)
        result.method = "unitigs"
        timings["unitig_compaction"] = time.perf_counter() - t0
        warnings.append(f"no Eulerian path ({euler.reason}); reporting {len(result.contigs)} unitigs")
    else:
        warnings.append("no k-mers survived (k larger than every read, or coverage filter too strict)")
    return result


def k_sweep(reads, ks, min_coverage: int = 1, clean: bool = False) -> list:
    """Assemble the same reads with several k values (the 'choosing k' experiment)."""
    rows = []
    for k in ks:
        r = assemble(reads, k, min_coverage=min_coverage, clean=clean)
        rows.append({
            "k": k,
            "vertices": r.graph.num_vertices(),
            "edges": r.graph.num_edges(),
            "method": r.method,
            "contigs": len(r.contigs),
            "longest": len(r.longest_contig),
            "ambiguous": r.ambiguous,
            "sequence": r.sequence,
        })
    return rows
