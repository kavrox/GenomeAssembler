"""Command-line interface for the In-Memory Genomic Sequence Assembler & Mutation Analyzer.

Every --reads / --ref / --query value may be a FASTA / FASTQ / plain-text file
OR a literal DNA string, e.g.

    python main.py pipeline --reads ATGC GCGAT GATAC -k 3 --ref ATGCGATAC --query ATGCAATAC
    python main.py assemble --reads data/error_reads.fasta -k 3 --min-coverage 2
    python main.py compare  --ref data/reference.fasta --query data/query.fasta
    python main.py ksweep   --reads data/reads.fasta --k 2-5
    python main.py search   --text ATGCGATAC ATGCAATAC --pattern ATAC
    python main.py tree     BANANA
    python main.py simulate --length 2000 --snv 3 --ins 1 --del 1 --out data/sim
    python main.py benchmark
    python main.py demo [k]
    python main.py gui
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import tracemalloc

from genome import report as rpt
from genome.assembler import assemble, k_sweep
from genome.io import format_fasta, load_reads, load_sequence, read_records, write_text
from genome.mutations import MutationAnalyzer
from genome.suffix_tree import GeneralizedSuffixTree


# ---------------------------------------------------------------- helpers
def parse_k_values(text: str) -> list:
    """'2-6' -> [2..6], '3,5,7' -> [3,5,7], '3,5-7' -> [3,5,6,7]."""
    ks = []
    for part in text.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = (int(x) for x in part.split("-", 1))
            if lo > hi:
                raise ValueError(f"bad k range {part!r}")
            ks.extend(range(lo, hi + 1))
        elif part:
            ks.append(int(part))
    if not ks or min(ks) < 2:
        raise ValueError("k values must be integers >= 2")
    return ks


def _load_reads(args):
    records, dropped = load_reads(args.reads, getattr(args, "min_quality", 0) or 0)
    if not records:
        raise ValueError("no reads left after quality filtering")
    return records, dropped


def _emit(text: str, collected: list):
    print(text)
    collected.append(text)


def _save_outputs(out_dir, files: dict):
    if not out_dir:
        return
    for name, text in files.items():
        write_text(os.path.join(out_dir, name), text)
    print(f"\nwrote {', '.join(sorted(files))} to {os.path.abspath(out_dir)}")


def _assembly_kwargs(args) -> dict:
    return {"min_coverage": args.min_coverage, "clean": args.clean, "multigraph": args.multigraph}


# ---------------------------------------------------------------- commands
def cmd_assemble(args) -> int:
    records, dropped = _load_reads(args)
    out = []
    if dropped:
        _emit(f"dropped {len(dropped)} read(s) with mean quality < {args.min_quality}", out)
    result = assemble(records, args.k, **_assembly_kwargs(args))
    _emit(rpt.format_assembly(result, verbose=not args.quiet), out)
    if args.profile:
        _emit("\n" + rpt.format_timings(result.timings), out)
    contigs = [(f"contig{i}_len{len(c)}", c) for i, c in enumerate(result.contigs, start=1)]
    _save_outputs(args.out, {"contigs.fasta": format_fasta(contigs) if contigs else "",
                             "graph.dot": result.graph.to_dot(result.path),
                             "assembly_report.txt": "\n".join(out) + "\n"})
    return 0 if result.contigs else 1


def cmd_compare(args) -> int:
    ref = load_sequence(args.ref, "reference")
    qry = load_sequence(args.query, "query")
    analyzer = MutationAnalyzer(ref.sequence, qry.sequence, args.min_anchor, args.algorithm)
    report = analyzer.report()
    out = []
    _emit(rpt.format_mutation_report(report, show_alignment=not args.no_alignment), out)
    if args.profile:
        _emit("\n" + rpt.format_timings(report.timings), out)
    _save_outputs(args.out, {"variants.tsv": rpt.variants_tsv(report.variants),
                             "mutation_report.txt": "\n".join(out) + "\n"})
    return 0


def cmd_pipeline(args) -> int:
    if args.profile:
        tracemalloc.start()
    t_start = time.perf_counter()
    records, dropped = _load_reads(args)
    reads = [r.sequence for r in records]
    ref = load_sequence(args.ref, "reference") if args.ref else None
    qry = load_sequence(args.query, "query") if args.query else None
    out = []

    _emit(rpt.banner("STEP 1  INPUT"), out)
    _emit(f"reads     : {len(reads)} (lengths {min(map(len, reads))}-{max(map(len, reads))} bp)"
          + (f", {len(dropped)} dropped for low quality" if dropped else ""), out)
    _emit(f"k         : {args.k}", out)
    _emit(f"reference : {rpt.shorten(ref.sequence) if ref else '-'}", out)
    _emit(f"query     : {rpt.shorten(qry.sequence) if qry else '(assembled sequence)'}", out)
    if len(reads) <= 20:
        _emit(rpt.banner("STEP 2  k-MERS  (sliding window, L-k+1 per read)"), out)
        _emit(rpt.format_kmers(reads, args.k), out)

    result = assemble(reads, args.k, **_assembly_kwargs(args))
    if len(result.freq) <= 60:
        _emit(f"\nSTEP 3  FREQUENCY TABLE: {rpt.format_frequency(result.freq)}", out)
    _emit(rpt.format_assembly(result, verbose=not args.quiet), out)

    sample = qry.sequence if qry else None
    if ref and result.contigs:
        if result.success:
            _emit(f"\nassembled sequence matches reference: {result.sequence == ref.sequence}", out)
            if sample is None:
                sample = result.sequence
        else:
            _emit("\nassembly is fragmented (no single Eulerian path): "
                  + ("comparing the given --query instead" if sample else
                     "mutation analysis skipped - fix k / --min-coverage / --clean, or give --query"), out)

    report = None
    if ref and sample:
        analyzer = MutationAnalyzer(ref.sequence, sample, args.min_anchor, args.algorithm)
        report = analyzer.report()
        _emit(rpt.format_mutation_report(report, show_alignment=not args.no_alignment), out)
    elif not ref:
        _emit("\n(no --ref given: mutation analysis skipped)", out)

    if args.profile:
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        _emit(rpt.banner("PROFILE"), out)
        _emit(rpt.format_timings(result.timings, report.timings if report else {}), out)
        _emit(f"total wall time: {rpt.format_seconds(time.perf_counter() - t_start)}   "
              f"peak traced memory: {peak / 1024:.1f} KiB", out)

    files = {"pipeline_report.txt": "\n".join(out) + "\n",
             "graph.dot": result.graph.to_dot(result.path)}
    if result.contigs:
        files["contigs.fasta"] = format_fasta(
            [(f"contig{i}_len{len(c)}", c) for i, c in enumerate(result.contigs, start=1)])
    if report:
        files["variants.tsv"] = rpt.variants_tsv(report.variants)
    _save_outputs(args.out, files)
    return 0


def cmd_ksweep(args) -> int:
    records, _ = _load_reads(args)
    ks = parse_k_values(args.k)
    rows = k_sweep([r.sequence for r in records], ks, args.min_coverage, args.clean)
    print(rpt.banner("k SWEEP - same reads, different k"))
    print(rpt.format_ksweep(rows))
    return 0


def _texts(tokens) -> list:
    """Files -> their records; anything else -> literal string (any alphabet)."""
    out = []
    for tok in tokens:
        if os.path.isfile(tok):
            out.extend(r.sequence for r in read_records(tok))
        else:
            out.append(tok)
    return out


def cmd_search(args) -> int:
    texts = _texts(args.text)
    tree = GeneralizedSuffixTree(texts, args.algorithm)
    print(f"GST over {len(texts)} sequence(s): {len(tree.text)} characters, {tree.count_nodes()} nodes")
    for p in args.pattern or []:
        occ = tree.occurrences(p)
        shown = ", ".join(f"seq{s + 1}:{o + 1}" for s, o in occ) or "not found"
        print(f"pattern {p!r}: {len(occ)} occurrence(s) (1-based)  {shown}")
    if len(texts) >= 2:
        lcs, pos = tree.lcs_with_positions()
        where = ", ".join(f"seq{s + 1}:{o + 1}" for s, o in pos.items())
        print(f"longest common substring: {rpt.shorten(lcs, 60) or '-'} (len {len(lcs)}) at {where or '-'}")
    for i in range(len(texts)):
        lrs = tree.longest_repeated_substring(i)
        print(f"longest repeated substring in seq{i + 1}: {rpt.shorten(lrs, 60) or '-'} (len {len(lrs)})")
    return 0


def cmd_tree(args) -> int:
    texts = _texts(args.text)
    tree = GeneralizedSuffixTree(texts, args.algorithm)
    joined = "".join(s + tree.sentinels[i] for i, s in enumerate(texts))
    print(f"suffix tree of {rpt.shorten(joined, 60)}  ({tree.count_nodes()} nodes, "
          f"{tree.count_leaves()} leaves, {args.algorithm})")
    print(tree.dump(args.max_lines))
    return 0


def cmd_simulate(args) -> int:
    from genome.simulate import mutate, random_genome, simulate_reads
    genome = random_genome(args.length, args.seed)
    reads = simulate_reads(genome, args.read_length, args.coverage, args.error_rate, args.seed)
    query, truth = mutate(genome, args.snv, args.ins, args.dels, seed=args.seed)
    files = {
        "reference.fasta": format_fasta([("reference", genome)]),
        "query.fasta": format_fasta([("query_mutated", query)]),
        "reads.fasta": format_fasta([(f"read{i}", r) for i, r in enumerate(reads, start=1)]),
        "truth.tsv": "type\tposition\tref\talt\n" + "".join(
            f"{t}\t{p}\t{r or '-'}\t{a or '-'}\n" for t, p, r, a in truth),
    }
    print(f"genome {args.length} bp, {len(reads)} reads of {args.read_length} bp "
          f"(~{args.coverage}x, error rate {args.error_rate}), {len(truth)} planted mutations:")
    for t in truth:
        print("  ", t)
    _save_outputs(args.out, files)
    return 0


def cmd_benchmark(args) -> int:
    from genome import benchmark as bm
    quick = args.quick
    g = lambda r, key: "-" if r.get(key + "_growth") is None else f"{r[key + '_growth']:.1f}x"
    ms = lambda t: f"{t * 1000:.1f}"

    print(rpt.banner("1. ASSEMBLY  (k=31, 100 bp reads, 10x coverage)  - expect ~2x per doubling"))
    rows = bm.assembly_scaling((500, 1000, 2000) if quick else (1000, 2000, 4000, 8000))
    print(rpt.format_table(
        ["genome n", "K kmers", "E", "count ms", "growth", "build ms", "growth", "euler ms", "growth", "correct"],
        [(r["n"], r["K"], r["E"], ms(r["count"]), g(r, "count"), ms(r["build"]), g(r, "build"),
          ms(r["euler"]), g(r, "euler"), r["correct"]) for r in rows]))

    print(rpt.banner("2. GST BUILD, random DNA  (two strings of length n)"))
    rows = bm.gst_scaling((500, 1000, 2000) if quick else (1000, 2000, 4000, 8000))
    print(rpt.format_table(["n", "text", "nodes", "naive ms", "growth", "ukkonen ms", "growth"],
                           [(r["n"], r["text"], r["nodes"], ms(r["naive"]), g(r, "naive"),
                             ms(r["ukkonen"]), g(r, "ukkonen")) for r in rows]))

    print(rpt.banner("3. GST BUILD, worst case A^n  - naive O(n^2) ~4x vs Ukkonen O(n) ~2x"))
    rows = bm.gst_scaling((100, 200, 400) if quick else (250, 500, 1000, 2000), worst_case=True)
    print(rpt.format_table(["n", "text", "nodes", "naive ms", "growth", "ukkonen ms", "growth"],
                           [(r["n"], r["text"], r["nodes"], ms(r["naive"]), g(r, "naive"),
                             ms(r["ukkonen"]), g(r, "ukkonen")) for r in rows]))

    print(rpt.banner("4. GST SEARCH  (pattern length m=20)  - expect flat: O(m), not O(n)"))
    rows = bm.search_scaling((1000, 4000) if quick else (1000, 4000, 16000))
    print(rpt.format_table(["text n", "m", "us per query"],
                           [(r["n"], r["m"], f"{r['per_query_us']:.2f}") for r in rows]))
    return 0


def cmd_demo(args) -> int:
    import demo
    demo.run(args.k)
    return 0


def cmd_gui(args) -> int:
    import gui
    gui.main()
    return 0


# ---------------------------------------------------------------- argument parser
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="main.py",
        description="In-Memory Genomic Sequence Assembler & Mutation Analyzer "
                    "(De Bruijn graph + Generalized Suffix Tree)")
    sub = p.add_subparsers(dest="command", metavar="command")

    def reads_opts(sp, k_required=True):
        sp.add_argument("--reads", nargs="+", required=True,
                        help="read files (FASTA/FASTQ/plain) and/or literal DNA reads")
        if k_required:
            sp.add_argument("-k", type=int, required=True, help="k-mer length (>= 2)")
        sp.add_argument("--min-coverage", type=int, default=1,
                        help="drop k-mers seen fewer times (error filter, default 1)")
        sp.add_argument("--clean", action="store_true", help="also remove tips and pop bubbles")
        sp.add_argument("--min-quality", type=float, default=0,
                        help="drop FASTQ reads with lower mean Phred quality")
        sp.add_argument("--multigraph", action="store_true",
                        help="one edge per k-mer occurrence (single known sequence only)")

    def compare_opts(sp):
        sp.add_argument("--min-anchor", type=int, default=3,
                        help="shortest shared block used as an anchor (default 3)")
        sp.add_argument("--algorithm", choices=["naive", "ukkonen"], default="naive",
                        help="suffix tree construction algorithm")
        sp.add_argument("--no-alignment", action="store_true", help="hide the alignment view")

    sp = sub.add_parser("assemble", help="assemble reads into contig(s)")
    reads_opts(sp)
    sp.add_argument("--quiet", action="store_true", help="hide the graph listing")
    sp.add_argument("--profile", action="store_true", help="show stage timings")
    sp.add_argument("--out", help="directory for contigs.fasta, graph.dot, report")
    sp.set_defaults(func=cmd_assemble)

    sp = sub.add_parser("compare", help="mutation analysis of reference vs query")
    sp.add_argument("--ref", required=True, help="reference file or DNA string")
    sp.add_argument("--query", required=True, help="query file or DNA string")
    compare_opts(sp)
    sp.add_argument("--profile", action="store_true", help="show stage timings")
    sp.add_argument("--out", help="directory for variants.tsv and the report")
    sp.set_defaults(func=cmd_compare)

    sp = sub.add_parser("pipeline", help="assemble reads, then compare with a reference")
    reads_opts(sp)
    sp.add_argument("--ref", help="reference file or DNA string")
    sp.add_argument("--query", help="query (default: the assembled sequence)")
    compare_opts(sp)
    sp.add_argument("--quiet", action="store_true", help="hide the graph listing")
    sp.add_argument("--profile", action="store_true", help="show timings and peak memory")
    sp.add_argument("--out", help="directory for all output files")
    sp.set_defaults(func=cmd_pipeline)

    sp = sub.add_parser("ksweep", help="assemble the same reads with several k values")
    reads_opts(sp, k_required=False)
    sp.add_argument("--k", required=True, help="k values, e.g. 2-6 or 3,5,7")
    sp.set_defaults(func=cmd_ksweep)

    sp = sub.add_parser("search", help="pattern search / LCS / repeats with a suffix tree")
    sp.add_argument("--text", nargs="+", required=True, help="files or literal strings")
    sp.add_argument("--pattern", nargs="*", help="patterns to locate")
    sp.add_argument("--algorithm", choices=["naive", "ukkonen"], default="naive")
    sp.set_defaults(func=cmd_search)

    sp = sub.add_parser("tree", help="print the (generalized) suffix tree")
    sp.add_argument("text", nargs="+", help="files or literal strings, e.g. BANANA")
    sp.add_argument("--algorithm", choices=["naive", "ukkonen"], default="naive")
    sp.add_argument("--max-lines", type=int, default=400)
    sp.set_defaults(func=cmd_tree)

    sp = sub.add_parser("simulate", help="generate a random genome, reads and a mutated query")
    sp.add_argument("--length", type=int, default=1000)
    sp.add_argument("--read-length", type=int, default=60)
    sp.add_argument("--coverage", type=float, default=10)
    sp.add_argument("--error-rate", type=float, default=0.0)
    sp.add_argument("--snv", type=int, default=2)
    sp.add_argument("--ins", type=int, default=1)
    sp.add_argument("--del", dest="dels", type=int, default=1)
    sp.add_argument("--seed", type=int, default=42)
    sp.add_argument("--out", default="data/sim")
    sp.set_defaults(func=cmd_simulate)

    sp = sub.add_parser("benchmark", help="measure how running time grows with input size")
    sp.add_argument("--quick", action="store_true", help="smaller sizes")
    sp.set_defaults(func=cmd_benchmark)

    sp = sub.add_parser("demo", help="the 8-step walk-through on ATGCGATAC")
    sp.add_argument("k", nargs="?", type=int, default=3)
    sp.set_defaults(func=cmd_demo)

    sp = sub.add_parser("gui", help="open the graphical interface")
    sp.set_defaults(func=cmd_gui)
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    try:
        return args.func(args)
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    sys.exit(main())
