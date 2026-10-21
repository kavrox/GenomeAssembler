"""Part 3 - Empirical complexity check.

Doubles the input size and measures how the running time of each stage grows:
    ~2x per doubling  -> linear     O(n)
    ~4x per doubling  -> quadratic  O(n^2)
    ~1x               -> independent of n (e.g. suffix-tree search, O(m))
"""

from __future__ import annotations

import random
import time

from .assembler import assemble
from .simulate import random_genome, simulate_reads
from .suffix_tree import GeneralizedSuffixTree


def _time(fn, repeat: int = 1):
    best = None
    for _ in range(repeat):
        t0 = time.perf_counter()
        value = fn()
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
    return best, value


def _ratio(rows, key):
    for prev, cur in zip(rows, rows[1:]):
        cur[key + "_growth"] = cur[key] / prev[key] if prev[key] > 0 else None
    if rows:
        rows[0][key + "_growth"] = None
    return rows


def assembly_scaling(sizes=(1000, 2000, 4000, 8000), k: int = 31, read_length: int = 100,
                     coverage: float = 10, seed: int = 7) -> list:
    rows = []
    for n in sizes:
        genome = random_genome(n, seed)
        reads = simulate_reads(genome, read_length, coverage, seed=seed)
        result = assemble(reads, k, check_ambiguity=False)
        t = result.timings
        rows.append({
            "n": n, "K": result.cleaning["total_kmers"], "E": result.graph.num_edges(),
            "count": t["kmer_counting"], "build": t["graph_construction"],
            "euler": t["eulerian_traversal"], "correct": result.sequence == genome,
        })
    for key in ("count", "build", "euler"):
        _ratio(rows, key)
    return rows


def gst_scaling(sizes=(1000, 2000, 4000, 8000), worst_case: bool = False, seed: int = 7) -> list:
    """Build time of the GST over two strings of length n each."""
    rows = []
    for n in sizes:
        if worst_case:
            seqs = ["A" * n, "A" * n]
        else:
            seqs = [random_genome(n, seed), random_genome(n, seed + 1)]
        naive, _ = _time(lambda: GeneralizedSuffixTree(seqs, "naive"))
        ukk, tree = _time(lambda: GeneralizedSuffixTree(seqs, "ukkonen"))
        rows.append({"n": n, "text": len(tree.text), "nodes": tree.count_nodes(),
                     "naive": naive, "ukkonen": ukk})
    for key in ("naive", "ukkonen"):
        _ratio(rows, key)
    return rows


def search_scaling(sizes=(1000, 4000, 16000), pattern_length: int = 20, queries: int = 2000,
                   seed: int = 7) -> list:
    """Pattern search cost should depend on the pattern length m, not the text size."""
    rng = random.Random(seed)
    rows = []
    for n in sizes:
        text = random_genome(n, seed)
        tree = GeneralizedSuffixTree([text], "ukkonen")
        patterns = []
        for _ in range(queries):
            s = rng.randint(0, n - pattern_length)
            patterns.append(text[s:s + pattern_length])
        t, _ = _time(lambda: [tree.contains(p) for p in patterns])
        rows.append({"n": n, "m": pattern_length, "per_query_us": t / queries * 1e6})
    return rows
