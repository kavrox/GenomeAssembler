"""Part 1 - Report generator: turns results into readable text / TSV.

Presentation only - no algorithmic logic lives here. Output is plain ASCII so
it prints correctly on any Windows console.
"""

from __future__ import annotations

from .kmers import generate_kmers, kmer_count_for_length

RULE = "=" * 64
NOTE = "NOTE: candidate sequence differences only - not a clinical interpretation."


def banner(title: str) -> str:
    return f"\n{RULE}\n{title}\n{RULE}"


def format_table(headers, rows) -> str:
    cells = [[str(h) for h in headers]] + [["" if c is None else str(c) for c in r] for r in rows]
    widths = [max(len(r[i]) for r in cells) for i in range(len(headers))]
    line = lambda r: "  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip()
    out = [line(cells[0]), "  ".join("-" * w for w in widths)]
    out.extend(line(r) for r in cells[1:])
    return "\n".join(out)


def format_seconds(t: float) -> str:
    return f"{t * 1000:.2f} ms" if t < 1 else f"{t:.2f} s"


def shorten(seq: str, width: int = 70) -> str:
    if seq is None:
        return "-"
    return seq if len(seq) <= width else f"{seq[:width // 2]}...{seq[-(width // 2 - 3):]} ({len(seq)} bp)"


# ---------------------------------------------------------------- assembly
def format_kmers(reads, k: int, max_reads: int = 20) -> str:
    lines = []
    for r in reads[:max_reads]:
        kmers = generate_kmers(r, k)
        shown = kmers if len(kmers) <= 12 else kmers[:12] + ["..."]
        lines.append(f"{shorten(r, 24):<24} L={len(r):<4} -> {kmer_count_for_length(len(r), k)} k-mers: {shown}")
    if len(reads) > max_reads:
        lines.append(f"... {len(reads) - max_reads} more reads")
    return "\n".join(lines)


def format_frequency(freq, max_items: int = 40) -> str:
    items = sorted(freq.items())
    shown = ", ".join(f"{k}:{c}" for k, c in items[:max_items])
    more = f", ... ({len(items)} distinct)" if len(items) > max_items else ""
    return "{" + shown + more + "}"


def format_graph(graph, max_vertices: int = 40) -> str:
    lines = []
    for v in sorted(graph.vertices)[:max_vertices]:
        targets = ", ".join(f"{e.dst}[{e.kmer} x{e.coverage}]" for e in graph.out[v]) or "-"
        mark = "   <- branch" if graph.out_degree(v) > 1 else ""
        lines.append(f"{v:<{graph.k + 1}} in={graph.in_degree(v)} out={graph.out_degree(v)}  ->  {targets}{mark}")
    if graph.num_vertices() > max_vertices:
        lines.append(f"... {graph.num_vertices() - max_vertices} more vertices")
    s = graph.summary()
    lines.append(f"vertices={s['vertices']} edges={s['edges']} "
                 f"branching={s['branching_vertices']} components={s['components']}")
    return "\n".join(lines)


def format_path(path, max_vertices: int = 16) -> str:
    if not path:
        return "-"
    if len(path) <= max_vertices:
        return " -> ".join(path)
    return (" -> ".join(path[:6]) + f" -> ... ({len(path) - 9} more) ... -> "
            + " -> ".join(path[-3:]))


def format_assembly(result, verbose: bool = True) -> str:
    out = [banner(f"ASSEMBLY  (k={result.k}, reads={len(result.reads)})")]
    c = result.cleaning
    out.append(f"k-mers: {c['total_kmers']} total, {c['distinct_kmers']} distinct")
    if c["low_coverage_removed"] or c["tip_edges_removed"] or c["bubble_edges_removed"]:
        out.append(f"cleaning: {c['low_coverage_removed']} low-coverage k-mers, "
                   f"{c['tip_edges_removed']} tip edges, {c['bubble_edges_removed']} bubble edges removed")
    if verbose:
        out.append("\nDe Bruijn graph (vertices=(k-1)-mers, edges=k-mers):")
        out.append(format_graph(result.graph))
    out.append(f"\nEuler check: {'OK' if result.euler.ok else 'FAILED'} - {result.euler.reason}")
    if result.path:
        out.append(f"Eulerian {('cycle' if result.euler.kind == 'cycle' else 'path')}: {format_path(result.path)}")
    out.append(f"\nmethod : {result.method}")
    if result.success:
        seq = result.sequence
        shown = seq if len(seq) <= 200 else shorten(seq, 100) + "  (full sequence: use --out)"
        out.append(f"assembled sequence ({len(seq)} bp): {shown}")
        if result.ambiguous is True:
            out.append(f"AMBIGUOUS - other sequences with the same k-mers: {result.alternatives}")
        elif result.ambiguous is False:
            out.append("unique: no other Euler path spells a different sequence")
    else:
        out.append(f"contigs: {len(result.contigs)}")
        for i, ctg in enumerate(result.contigs, start=1):
            out.append(f"  contig{i} ({len(ctg)} bp): {shorten(ctg)}")
    for w in result.warnings:
        out.append(f"warning: {w}")
    return "\n".join(out)


def format_ksweep(rows) -> str:
    table = [(r["k"], r["vertices"], r["edges"], r["method"], r["contigs"], r["longest"],
              {True: "yes", False: "no", None: "-"}[r["ambiguous"]], shorten(r["sequence"], 40))
             for r in rows]
    return format_table(["k", "V", "E", "result", "contigs", "longest", "ambiguous", "sequence"], table)


# ---------------------------------------------------------------- mutations
def format_alignment(report, width: int = 60) -> str:
    top, mid, bot = report.alignment()
    out = []
    for i in range(0, len(top), width):
        out.append(f"REF   {top[i:i + width]}")
        out.append(f"      {mid[i:i + width]}")
        out.append(f"QUERY {bot[i:i + width]}")
        out.append("")
    return "\n".join(out).rstrip()


def format_variants(variants) -> str:
    if not variants:
        return "no differences: reference and query are identical"
    rows = [(i, v.type, v.position, v.ref or "-", v.alt or "-", v.hgvs, v.context)
            for i, v in enumerate(variants, start=1)]
    return format_table(["#", "type", "pos", "ref", "alt", "notation", "context"], rows)


def format_mutation_report(report, show_alignment: bool = True) -> str:
    out = [banner(f"MUTATION ANALYSIS  (GST, {report.algorithm} construction)")]
    out.append(f"reference: {len(report.reference)} bp   query: {len(report.query)} bp")
    out.append(f"GST over ref#query$: {report.gst_text_length} characters, {report.gst_nodes} nodes")
    if report.lcs:
        p = report.lcs_positions
        out.append(f"longest common substring: {shorten(report.lcs, 50)} "
                   f"(len {len(report.lcs)}, ref pos {p[0] + 1}, query pos {p[1] + 1})")
    out.append(f"anchors: {len(report.anchors)}   identity (anchored fraction): {report.identity:.1%}")
    if report.low_similarity:
        out.append("warning: LOW SIMILARITY - sequences share little; the variant list is not meaningful")
    if show_alignment and len(report.reference) + len(report.query) <= 4000:
        out.append("\n" + format_alignment(report))
    out.append("")
    out.append(format_variants(report.variants))
    counts = report.counts
    if counts:
        out.append("\nsummary: " + ", ".join(f"{n} {t.lower()}" for t, n in sorted(counts.items())))
    out.append(NOTE)
    return "\n".join(out)


def variants_tsv(variants) -> str:
    lines = ["type\tposition\tref\talt\thgvs\tcontext"]
    lines.extend(f"{v.type}\t{v.position}\t{v.ref or '-'}\t{v.alt or '-'}\t{v.hgvs}\t{v.context}"
                 for v in variants)
    return "\n".join(lines) + "\n"


def format_timings(*timing_dicts) -> str:
    rows = []
    for d in timing_dicts:
        rows.extend((name.replace("_", " "), format_seconds(t)) for name, t in d.items())
    return format_table(["stage", "time"], rows)
