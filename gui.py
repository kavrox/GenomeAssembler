"""Part 2 - Tkinter GUI: draws the De Bruijn graph, the Euler path, the alignment
and the variant table. Standard library only.

    python gui.py
"""

from __future__ import annotations

import math
import random
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from genome import report as rpt
from genome.assembler import assemble, k_sweep
from genome.io import parse_text
from genome.mutations import MutationAnalyzer

EXAMPLES = {
    "Textbook: 3 reads, SNV (k=3)": dict(
        reads="ATGC\nGCGAT\nGATAC", k=3, ref="ATGCGATAC", query="ATGCAATAC", mincov=1, clean=False),
    "Sequencing error -> bubble (k=3)": dict(
        reads="ATGCGATAC\nATGCGATAC\nATGCGATAC\nATGCGTTAC", k=3, ref="ATGCGATAC", query="",
        mincov=1, clean=False),
    "Error removed by cleaning (k=3)": dict(
        reads="ATGCGATAC\nATGCGATAC\nATGCGATAC\nATGCGTTAC", k=3, ref="ATGCGATAC", query="",
        mincov=1, clean=True),
    "k too large -> fragments (k=4)": dict(
        reads="ATGC\nGCGAT\nGATAC", k=4, ref="ATGCGATAC", query="", mincov=1, clean=False),
    "Repeat ambiguity TTGTGAGTA (k=3)": dict(
        reads="TTGTGAG\nGTGAGTA", k=3, ref="TTGTGAGTA", query="", mincov=1, clean=False),
    "SNV + insertion + deletion (k=5)": dict(
        reads="ATGCGATACCGTAGGCT\nCCGTAGGCTTACGATCA", k=5,
        ref="ATGCGATACCGTAGGCTTACGATCA", query="ATGCAATACCGGTAGGCTTACGTCA", mincov=1, clean=False),
}

COL = {
    "bg": "#f7f8fa", "node": "#dbe7ff", "node_line": "#1e2a44", "branch": "#fff1bf",
    "start": "#c9f2d9", "end": "#ffd6d6", "edge": "#7a8699", "path": "#1f6feb",
    "label": "#8a3b12", "match": "#1a7f37", "mismatch": "#cf222e", "indel": "#9a6700",
}


# ---------------------------------------------------------------- graph layout
def force_layout(graph, width, height, path=None, iterations=300, seed=1):
    """Fruchterman-Reingold spring layout (deterministic). Vertices on an Euler path
    start on a circle in walk order, which keeps the drawing readable."""
    verts = sorted(graph.vertices)
    n = len(verts)
    if n == 0:
        return {}
    order = []
    for v in (path or []) + verts:
        if v not in order:
            order.append(v)
    pos = {}
    for i, v in enumerate(order):
        a = 2 * math.pi * i / n
        pos[v] = [width / 2 + 0.35 * width * math.cos(a), height / 2 + 0.35 * height * math.sin(a)]
    if n == 1:
        return {verts[0]: (width / 2, height / 2)}
    rng = random.Random(seed)
    area = width * height
    kk = math.sqrt(area / n) * 0.75
    edges = [(e.src, e.dst) for e in graph.edges() if e.src != e.dst]
    temp = width / 8
    for _ in range(iterations):
        disp = {v: [0.0, 0.0] for v in verts}
        for i, a in enumerate(verts):
            for b in verts[i + 1:]:
                dx, dy = pos[a][0] - pos[b][0], pos[a][1] - pos[b][1]
                d = math.hypot(dx, dy) or rng.uniform(0.01, 0.1)
                f = kk * kk / d
                disp[a][0] += dx / d * f
                disp[a][1] += dy / d * f
                disp[b][0] -= dx / d * f
                disp[b][1] -= dy / d * f
        for a, b in edges:
            dx, dy = pos[a][0] - pos[b][0], pos[a][1] - pos[b][1]
            d = math.hypot(dx, dy) or 0.01
            f = d * d / kk
            disp[a][0] -= dx / d * f
            disp[a][1] -= dy / d * f
            disp[b][0] += dx / d * f
            disp[b][1] += dy / d * f
        for v in verts:
            dx, dy = disp[v]
            d = math.hypot(dx, dy) or 1
            pos[v][0] = min(width - 50, max(50, pos[v][0] + dx / d * min(d, temp)))
            pos[v][1] = min(height - 40, max(40, pos[v][1] + dy / d * min(d, temp)))
        temp *= 0.97
    return {v: tuple(p) for v, p in pos.items()}


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("In-Memory Genomic Sequence Assembler & Mutation Analyzer")
        self.geometry("1280x800")
        self.minsize(1000, 640)
        self.result = None
        self.mreport = None
        self._build_inputs()
        self._build_outputs()
        self.load_example(next(iter(EXAMPLES)))

    # ------------------------------------------------------------ layout
    def _build_inputs(self):
        left = ttk.Frame(self, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(left, text="Example").pack(anchor="w")
        self.example = ttk.Combobox(left, values=list(EXAMPLES), state="readonly", width=36)
        self.example.pack(fill=tk.X)
        self.example.bind("<<ComboboxSelected>>", lambda e: self.load_example(self.example.get()))

        head = ttk.Frame(left)
        head.pack(fill=tk.X, pady=(10, 0))
        ttk.Label(head, text="Reads (one per line, or FASTA/FASTQ)").pack(side=tk.LEFT)
        ttk.Button(head, text="Open...", command=self.open_reads).pack(side=tk.RIGHT)
        self.reads_text = tk.Text(left, width=40, height=12, font=("Consolas", 10))
        self.reads_text.pack(fill=tk.X)

        params = ttk.Frame(left)
        params.pack(fill=tk.X, pady=8)
        ttk.Label(params, text="k").grid(row=0, column=0, sticky="w")
        self.k_var = tk.IntVar(value=3)
        ttk.Spinbox(params, from_=2, to=99, textvariable=self.k_var, width=5).grid(row=0, column=1, sticky="w")
        ttk.Label(params, text="min coverage").grid(row=0, column=2, sticky="w", padx=(12, 0))
        self.mincov_var = tk.IntVar(value=1)
        ttk.Spinbox(params, from_=1, to=999, textvariable=self.mincov_var, width=5).grid(row=0, column=3, sticky="w")
        self.clean_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(params, text="remove tips / pop bubbles", variable=self.clean_var).grid(
            row=1, column=0, columnspan=4, sticky="w", pady=(4, 0))
        ttk.Label(params, text="suffix tree").grid(row=2, column=0, sticky="w", pady=(4, 0))
        self.algo_var = tk.StringVar(value="naive")
        ttk.Combobox(params, textvariable=self.algo_var, values=["naive", "ukkonen"],
                     state="readonly", width=9).grid(row=2, column=1, columnspan=2, sticky="w", pady=(4, 0))

        for label, attr, cmd in (("Reference", "ref_text", self.open_ref),
                                 ("Query (blank = assembled sequence)", "query_text", self.open_query)):
            row = ttk.Frame(left)
            row.pack(fill=tk.X, pady=(6, 0))
            ttk.Label(row, text=label).pack(side=tk.LEFT)
            ttk.Button(row, text="Open...", command=cmd).pack(side=tk.RIGHT)
            widget = tk.Text(left, width=40, height=3, font=("Consolas", 10))
            widget.pack(fill=tk.X)
            setattr(self, attr, widget)

        buttons = ttk.Frame(left)
        buttons.pack(fill=tk.X, pady=10)
        ttk.Button(buttons, text="Run pipeline", command=self.run).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(buttons, text="k sweep 2-8", command=self.run_ksweep).pack(side=tk.LEFT, expand=True, fill=tk.X)

        self.status = ttk.Label(left, text="", wraplength=300, foreground="#444")
        self.status.pack(fill=tk.X)

    def _build_outputs(self):
        self.tabs = ttk.Notebook(self)
        self.tabs.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(0, 8), pady=8)

        frame = ttk.Frame(self.tabs)
        self.canvas = tk.Canvas(frame, background=COL["bg"], highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self.canvas.bind("<Configure>", lambda e: self.draw_graph())
        legend = ("yellow = branch vertex   green = path start   red = path end   "
                  "blue numbers = Euler walk order   label = k-mer x coverage")
        ttk.Label(frame, text=legend, foreground="#555").pack(anchor="w")
        self.tabs.add(frame, text="De Bruijn graph")

        self.assembly_out = self._text_tab("Assembly")

        frame = ttk.Frame(self.tabs)
        self.align_out = tk.Text(frame, height=10, font=("Consolas", 11), wrap="none")
        self.align_out.pack(fill=tk.X)
        for tag in ("match", "mismatch", "indel"):
            self.align_out.tag_configure(tag, foreground=COL[tag])
        self.align_out.tag_configure("marker", foreground="#8c959f")
        cols = ("type", "pos", "ref", "alt", "notation", "context")
        self.var_table = ttk.Treeview(frame, columns=cols, show="headings", height=8)
        for c, w in zip(cols, (120, 60, 80, 80, 140, 240)):
            self.var_table.heading(c, text=c)
            self.var_table.column(c, width=w, anchor="w")
        self.var_table.pack(fill=tk.BOTH, expand=True, pady=6)
        self.mut_summary = tk.Text(frame, height=7, font=("Consolas", 10), wrap="word")
        self.mut_summary.pack(fill=tk.X)
        self.tabs.add(frame, text="Mutations")

        self.tree_out = self._text_tab("Suffix tree")
        self.ksweep_out = self._text_tab("k sweep")

    def _text_tab(self, title):
        frame = ttk.Frame(self.tabs)
        text = tk.Text(frame, font=("Consolas", 10), wrap="none")
        sb = ttk.Scrollbar(frame, command=text.yview)
        text.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        text.pack(fill=tk.BOTH, expand=True)
        self.tabs.add(frame, text=title)
        return text

    # ------------------------------------------------------------ helpers
    @staticmethod
    def _set(widget, text):
        widget.configure(state="normal")
        widget.delete("1.0", tk.END)
        widget.insert("1.0", text)

    @staticmethod
    def _get(widget):
        return widget.get("1.0", tk.END).strip()

    def load_example(self, name):
        ex = EXAMPLES[name]
        self.example.set(name)
        self._set(self.reads_text, ex["reads"])
        self._set(self.ref_text, ex["ref"])
        self._set(self.query_text, ex["query"])
        self.k_var.set(ex["k"])
        self.mincov_var.set(ex["mincov"])
        self.clean_var.set(ex["clean"])
        self.run()

    def _open_into(self, widget, first_only):
        path = filedialog.askopenfilename(filetypes=[("Sequence files", "*.fasta *.fa *.fastq *.fq *.txt"),
                                                     ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as fh:
                records = parse_text(fh.read())
        except (OSError, ValueError) as exc:
            messagebox.showerror("Could not read file", str(exc))
            return
        if first_only:
            self._set(widget, records[0].sequence if records else "")
        else:
            self._set(widget, "\n".join(r.sequence for r in records))

    def open_reads(self):
        self._open_into(self.reads_text, False)

    def open_ref(self):
        self._open_into(self.ref_text, True)

    def open_query(self):
        self._open_into(self.query_text, True)

    def _reads(self):
        text = self._get(self.reads_text)
        if not text:
            raise ValueError("enter at least one read")
        return [r.sequence for r in parse_text(text)]

    # ------------------------------------------------------------ actions
    def run(self):
        try:
            reads = self._reads()
            k = int(self.k_var.get())
            result = assemble(reads, k, min_coverage=int(self.mincov_var.get()), clean=self.clean_var.get())
            ref_text = "".join(self._get(self.ref_text).split())
            query_text = "".join(self._get(self.query_text).split())
            mreport = None
            sample = query_text or result.sequence
            if ref_text and sample:
                mreport = MutationAnalyzer(ref_text, sample, algorithm=self.algo_var.get()).report()
        except (ValueError, tk.TclError) as exc:
            self.status.configure(text=f"Error: {exc}", foreground=COL["mismatch"])
            return
        self.result, self.mreport = result, mreport
        self._show_assembly(ref_text)
        self._show_mutations()
        self.draw_graph()
        msg = (f"Assembled {len(result.sequence)} bp via {result.method}" if result.success
               else f"No Eulerian path - {len(result.contigs)} unitig(s)")
        if mreport:
            msg += f"; {len(mreport.variants)} candidate variant(s)"
        self.status.configure(text=msg, foreground="#1a7f37" if result.success else COL["indel"])

    def run_ksweep(self):
        try:
            rows = k_sweep(self._reads(), range(2, 9), int(self.mincov_var.get()), self.clean_var.get())
        except (ValueError, tk.TclError) as exc:
            self.status.configure(text=f"Error: {exc}", foreground=COL["mismatch"])
            return
        self._set(self.ksweep_out, "Same reads, different k\n\n" + rpt.format_ksweep(rows))
        self.tabs.select(4)

    def _show_assembly(self, ref_text):
        text = rpt.format_assembly(self.result)
        if ref_text and self.result.contigs:
            assembled = self.result.sequence or self.result.longest_contig
            text += f"\n\nassembled sequence matches reference: {assembled == ref_text}"
        text += "\n\n" + rpt.format_timings(self.result.timings)
        self._set(self.assembly_out, text)

    def _show_mutations(self):
        for item in self.var_table.get_children():
            self.var_table.delete(item)
        rep = self.mreport
        if rep is None:
            if self._get(self.ref_text) and not self.result.success:
                msg = ("Assembly is fragmented (no single Eulerian path), so there is no single sequence "
                       "to compare.\nFix k / min coverage / cleaning, or type a query sequence.")
            else:
                msg = "Enter a reference to run mutation analysis."
            self._set(self.align_out, msg)
            self._set(self.mut_summary, "")
            self._set(self.tree_out, "")
            return
        top, mid, bot = rep.alignment()
        self._set(self.align_out, "")
        width = 80
        for i in range(0, len(top), width):
            for label, row in (("REF   ", top), ("      ", mid), ("QUERY ", bot)):
                self.align_out.insert(tk.END, label)
                for j, ch in enumerate(row[i:i + width]):
                    m = mid[i + j]
                    if row is mid:
                        tag = "mismatch" if m == "x" else "marker"
                    else:
                        tag = "match" if m == "|" else ("mismatch" if m == "x" else "indel")
                    self.align_out.insert(tk.END, ch, tag)
                self.align_out.insert(tk.END, "\n")
            self.align_out.insert(tk.END, "\n")
        for v in rep.variants:
            self.var_table.insert("", tk.END, values=(v.type, v.position, v.ref or "-", v.alt or "-",
                                                      v.hgvs, v.context))
        p = rep.lcs_positions
        lines = [f"GST over ref#query$: {rep.gst_text_length} characters, {rep.gst_nodes} nodes ({rep.algorithm})",
                 f"longest common substring: {rpt.shorten(rep.lcs, 60)} (len {len(rep.lcs)}"
                 + (f", ref pos {p[0] + 1}, query pos {p[1] + 1})" if rep.lcs else ")"),
                 f"anchors: {len(rep.anchors)}   identity: {rep.identity:.1%}"]
        if rep.low_similarity:
            lines.append("LOW SIMILARITY - the variant list is not meaningful")
        if not rep.variants:
            lines.append("no differences: reference and query are identical")
        lines.append(rpt.NOTE)
        self._set(self.mut_summary, "\n".join(lines))
        if rep.gst_text_length <= 200:
            from genome.suffix_tree import GeneralizedSuffixTree
            tree = GeneralizedSuffixTree([rep.reference, rep.query], rep.algorithm)
            self._set(self.tree_out, f"Generalized suffix tree of {tree.text}\n\n{tree.dump()}")
        else:
            self._set(self.tree_out, f"Tree has {rep.gst_nodes} nodes - too large to print here "
                                     "(use: python main.py tree <files>).")

    # ------------------------------------------------------------ drawing
    def draw_graph(self):
        c = self.canvas
        c.delete("all")
        if self.result is None:
            return
        g = self.result.graph
        w, h = max(c.winfo_width(), 400), max(c.winfo_height(), 300)
        if g.num_vertices() == 0:
            c.create_text(w / 2, h / 2, text="Empty graph: no k-mers (is k longer than the reads?)")
            return
        if g.num_vertices() > 60:
            c.create_text(w / 2, h / 2, text=f"Graph too large to draw ({g.num_vertices()} vertices, "
                                             f"{g.num_edges()} edges) - see the Assembly tab.")
            return
        path = self.result.path or []
        pos = force_layout(g, w, h, path)
        steps = {}
        for i, (a, b) in enumerate(zip(path, path[1:]), start=1):
            steps.setdefault((a, b), []).append(i)
        radius = max(16, 5 * (g.k - 1) + 6)
        pair_seen = {}
        for e in sorted(g.edges(), key=lambda e: (e.src, e.kmer)):
            (x1, y1), (x2, y2) = pos[e.src], pos[e.dst]
            on_path = (e.src, e.dst) in steps
            colour = COL["path"] if on_path else COL["edge"]
            if e.src == e.dst:                                      # self loop
                c.create_oval(x1 - radius, y1 - 2.6 * radius, x1 + radius, y1 - 0.6 * radius,
                              outline=colour, width=2)
                c.create_text(x1, y1 - 2.9 * radius, text=f"{e.kmer} x{e.coverage}",
                              fill=COL["label"], font=("Consolas", 9))
                continue
            n_same = pair_seen.get((e.src, e.dst), 0)
            pair_seen[(e.src, e.dst)] = n_same + 1
            reverse = g.has_edge(e.dst, e.src)
            dx, dy = x2 - x1, y2 - y1
            d = math.hypot(dx, dy) or 1
            ux, uy = dx / d, dy / d
            sx, sy = x1 + ux * radius, y1 + uy * radius
            ex, ey = x2 - ux * radius, y2 - uy * radius
            bend = (25 if reverse else 0) + 22 * n_same
            mx, my = (sx + ex) / 2 - uy * bend, (sy + ey) / 2 + ux * bend
            c.create_line(sx, sy, mx, my, ex, ey, smooth=True, arrow=tk.LAST, arrowshape=(11, 13, 5),
                          fill=colour, width=2.5 if on_path else 1.5)
            label = f"{e.kmer} x{e.coverage}"
            if on_path and steps[(e.src, e.dst)]:
                label += f"  #{steps[(e.src, e.dst)].pop(0)}"
            tid = c.create_text(mx - uy * 12, my + ux * 12, text=label,
                                fill=COL["path"] if on_path else COL["label"], font=("Consolas", 9))
            bx = c.bbox(tid)
            c.tag_lower(c.create_rectangle(bx[0] - 2, bx[1] - 1, bx[2] + 2, bx[3] + 1,
                                           fill=COL["bg"], outline=""), tid)
        start = path[0] if path else None
        end = path[-1] if path else None
        for v, (x, y) in pos.items():
            fill = COL["node"]
            if g.out_degree(v) > 1 or g.in_degree(v) > 1:
                fill = COL["branch"]
            if v == start:
                fill = COL["start"]
            elif v == end:
                fill = COL["end"]
            c.create_oval(x - radius, y - radius, x + radius, y + radius, fill=fill,
                          outline=COL["node_line"], width=2)
            c.create_text(x, y, text=v, font=("Consolas", 10, "bold"))
        c.create_text(10, 10, anchor="nw", font=("Segoe UI", 10),
                      text=f"k={g.k}   V={g.num_vertices()}   E={g.num_edges()}   {self.result.method}")


def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
