"""Live demo (under 5 minutes): the 8 steps of the project on one tiny genome.

    python demo.py        k = 3  -> reconstructs ATGCGATAC, finds G>A at position 5
    python demo.py 4      k = 4  -> "k too large": the graph fragments, no Euler path
    python demo.py 2      k = 2  -> "k too small": a valid Euler path spells the WRONG genome
"""

from __future__ import annotations

import sys

from genome.assembler import assemble
from genome.euler import check_euler, path_to_sequence
from genome.kmers import generate_kmers
from genome.mutations import MutationAnalyzer
from genome.report import NOTE, banner, format_alignment, format_graph, format_path

READS = ["ATGC", "GCGAT", "GATAC"]
REFERENCE = "ATGCGATAC"
QUERY = "ATGCAATAC"


def run(k: int = 3, reads=READS, reference: str = REFERENCE, query: str = QUERY) -> None:
    print(banner("STEP 1  INPUT"))
    print(f"reads : {reads}")
    print(f"k     : {k}")
    print(f"ref   : {reference}")
    print(f"query : {query}")

    print(banner("STEP 2  k-MERS  (window slides 1 step; L-k+1 per read)"))
    for r in reads:
        km = generate_kmers(r, k)
        print(f"{r:<8} L={len(r)} -> {len(km)} k-mers: {km}")

    result = assemble(reads, k)
    print(f"\nSTEP 3  FREQUENCY TABLE: {dict(sorted(result.freq.items()))}")

    print(banner("STEP 4  DE BRUIJN GRAPH  (vertices=(k-1)-mers, edges=k-mers)"))
    print(format_graph(result.graph))
    ec = check_euler(result.graph)
    print(f"Euler start={ec.start} end={ec.end} conditions-ok={ec.ok}  ({ec.reason})")

    print(banner("STEP 5  EULERIAN PATH  (Hierholzer)"))
    if result.path:
        print(format_path(result.path))
    else:
        print(f"NO EULERIAN PATH - {ec.reason}")
        print(f"fallback unitigs: {result.contigs}")
        if k > 3:
            print("k is too large for these reads: neighbouring reads overlap by fewer than k-1 letters.")
        return

    assembled = path_to_sequence(result.path)
    print(banner("STEP 6  ASSEMBLED SEQUENCE"))
    print(f"assembled        : {assembled}")
    print(f"matches reference: {assembled == reference}")
    if assembled != reference:
        print("k is too small: short k-mers create accidental overlaps, so a valid Euler path "
              "spells the wrong sequence.")
    if result.ambiguous:
        print(f"other sequences using the same k-mers: {result.alternatives}")

    analyzer = MutationAnalyzer(reference, query)
    rep = analyzer.report()
    print(banner("STEP 7  GENERALIZED SUFFIX TREE  (ref # query $)"))
    print(f"nodes in tree: {rep.gst_nodes}   text length: {rep.gst_text_length}")
    print(f"longest common substring: {rep.lcs}")
    print(f"occurrences of 'ATAC' (string_id, offset): {analyzer.gst.occurrences('ATAC')}")

    print(banner("STEP 8  CANDIDATE MUTATIONS"))
    print(format_alignment(rep))
    print()
    if not rep.variants:
        print("no differences")
    for v in rep.variants:
        print(f"{v.type:<13} position={v.position}  ref={v.ref or '-'}  query={v.alt or '-'}  "
              f"{v.hgvs:<12} context: {v.context}")
    print(f"\n{NOTE}")


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
