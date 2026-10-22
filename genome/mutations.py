"""Part 4 - Mutation analysis with the generalized suffix tree.

Anchor-gap algorithm (a longest-common-block strategy, like Ratcliff-Obershelp):
  1. build a GST over ref + '#' + query + '$'
  2. the longest common substring (deepest node with leaves from both strings)
     is an ANCHOR - a stretch we trust is identical
  3. repeat on the part left of the anchor and the part right of it
  4. when no anchor of at least ``min_anchor`` letters is left, the leftover
     pair (ref gap, query gap) is a candidate difference and is classified:
        equal length        -> SUBSTITUTION(s)
        empty ref gap       -> INSERTION
        empty query gap     -> DELETION
        otherwise           -> trim, then INS/DEL or COMPLEX (delins)

This reports candidate sequence differences only. It knows nothing about
biology: it cannot say whether a difference is real, common (SNP) or harmful.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from .sequence import validate_dna
from .suffix_tree import GeneralizedSuffixTree


@dataclass
class Variant:
    type: str          # SUBSTITUTION | INSERTION | DELETION | COMPLEX
    position: int      # 1-based reference position (INSERTION: inserted before it)
    ref: str           # reference allele ("" for an insertion)
    alt: str           # query allele ("" for a deletion)
    context: str = ""

    def as_tuple(self):
        return (self.type, self.position, self.ref, self.alt)

    @property
    def hgvs(self) -> str:
        """HGVS-style genomic notation, e.g. g.5G>A, g.5del, g.5_6insG."""
        p = self.position
        if self.type == "SUBSTITUTION":
            return f"g.{p}{self.ref}>{self.alt}"
        if self.type == "INSERTION":
            return f"g.{p - 1}_{p}ins{self.alt}"
        span = f"{p}" if len(self.ref) == 1 else f"{p}_{p + len(self.ref) - 1}"
        if self.type == "DELETION":
            return f"g.{span}del"
        return f"g.{span}delins{self.alt}"


@dataclass
class Anchor:
    ref_start: int     # 0-based
    query_start: int   # 0-based
    length: int

    @property
    def ref_end(self):
        return self.ref_start + self.length

    @property
    def query_end(self):
        return self.query_start + self.length


# ---------------------------------------------------------------- gap classification
def common_prefix_len(a: str, b: str) -> int:
    n = 0
    while n < min(len(a), len(b)) and a[n] == b[n]:
        n += 1
    return n


def common_suffix_len(a: str, b: str, limit: int = None) -> int:
    limit = min(len(a), len(b)) if limit is None else limit
    n = 0
    while n < limit and a[len(a) - 1 - n] == b[len(b) - 1 - n]:
        n += 1
    return n


def classify_gap(rgap: str, qgap: str, ref_pos: int) -> list:
    """Classify ONE unmatched gap. ref_pos = 1-based reference index of the gap start."""
    if rgap == qgap:
        return []
    if not rgap:
        return [Variant("INSERTION", ref_pos, "", qgap)]
    if not qgap:
        return [Variant("DELETION", ref_pos, rgap, "")]
    if len(rgap) == len(qgap):                          # per-base mismatches
        return [Variant("SUBSTITUTION", ref_pos + i, a, b)
                for i, (a, b) in enumerate(zip(rgap, qgap)) if a != b]
    p = common_prefix_len(rgap, qgap)                   # unequal: trim, then re-test
    s = common_suffix_len(rgap[p:], qgap[p:])
    r2, q2 = rgap[p:len(rgap) - s], qgap[p:len(qgap) - s]
    if not r2 or not q2:
        return classify_gap(r2, q2, ref_pos + p)
    return [Variant("COMPLEX", ref_pos + p, r2, q2)]


def classify_difference(ref: str, qry: str) -> list:
    """Single-edit shortcut: trim common prefix/suffix (valid for ONE localised edit)."""
    if ref == qry:
        return []
    p = common_prefix_len(ref, qry)
    s = common_suffix_len(ref[p:], qry[p:])
    return classify_gap(ref[p:len(ref) - s], qry[p:len(qry) - s], p + 1)


# ---------------------------------------------------------------- anchor search
def anchor_alignment(ref: str, qry: str, min_anchor: int = 3, algorithm: str = "naive"):
    """Split ref/query into ordered anchors (identical blocks) and gaps.

    Iterative (explicit stack) so very different sequences cannot overflow the
    recursion limit. Returns (anchors, gaps) with gaps as (r0, r1, q0, q1).
    """
    anchors, gaps = [], []
    stack = [(0, len(ref), 0, len(qry))]
    while stack:
        r0, r1, q0, q1 = stack.pop()
        rs, qs = ref[r0:r1], qry[q0:q1]
        if rs == qs:
            if rs:
                anchors.append(Anchor(r0, q0, len(rs)))
            continue
        lcs, pos = "", {}
        if rs and qs:
            lcs, pos = GeneralizedSuffixTree([rs, qs], algorithm).lcs_with_positions()
        if len(lcs) < min_anchor:
            gaps.append((r0, r1, q0, q1))
            continue
        ri, qi = r0 + pos[0], q0 + pos[1]
        anchors.append(Anchor(ri, qi, len(lcs)))
        stack.append((ri + len(lcs), r1, qi + len(lcs), q1))
        stack.append((r0, ri, q0, qi))
    anchors.sort(key=lambda a: (a.ref_start, a.query_start))
    gaps.sort(key=lambda g: (g[0], g[2]))
    return anchors, gaps


def _context(ref: str, v: Variant, flank: int) -> str:
    i = v.position - 1
    left = ref[max(0, i - flank):i]
    right = ref[i + len(v.ref):i + len(v.ref) + flank]
    return f"{left}[{v.ref or '-'}>{v.alt or '-'}]{right}"


def find_variants(ref: str, qry: str, min_anchor: int = 3, algorithm: str = "naive",
                  flank: int = 4) -> list:
    """Candidate variants between reference and query, ordered by position."""
    if ref == qry:
        return []
    _, gaps = anchor_alignment(ref, qry, min_anchor, algorithm)
    variants = []
    for r0, r1, q0, q1 in gaps:
        variants.extend(classify_gap(ref[r0:r1], qry[q0:q1], r0 + 1))
    variants.sort(key=lambda v: (v.position, v.type))
    for v in variants:
        v.context = _context(ref, v, flank)
    return variants


def build_alignment(ref: str, qry: str, anchors: list, gaps: list):
    """Gapped alignment rows (ref, markers, query) from anchors and gaps.

    markers: '|' match, 'x' substitution, ' ' insertion/deletion column.
    """
    items = [("A", a.ref_start, a.query_start, a) for a in anchors]
    items += [("G", g[0], g[2], g) for g in gaps]
    items.sort(key=lambda t: (t[1], t[2], t[0] == "A"))
    top, mid, bot = [], [], []
    for kind, _, _, item in items:
        if kind == "A":
            block = ref[item.ref_start:item.ref_end]
            top.append(block)
            mid.append("|" * len(block))
            bot.append(block)
            continue
        r, q = ref[item[0]:item[1]], qry[item[2]:item[3]]
        if len(r) == len(q):
            top.append(r)
            bot.append(q)
            mid.append("".join("|" if a == b else "x" for a, b in zip(r, q)))
            continue
        p = common_prefix_len(r, q)
        s = common_suffix_len(r[p:], q[p:])
        r_mid, q_mid = r[p:len(r) - s], q[p:len(q) - s]
        width = max(len(r_mid), len(q_mid))
        r_row = r[:p] + r_mid.ljust(width, "-") + r[len(r) - s:]
        q_row = q[:p] + q_mid.ljust(width, "-") + q[len(q) - s:]
        top.append(r_row)
        bot.append(q_row)
        mid.append("".join("|" if a == b else (" " if "-" in (a, b) else "x")
                           for a, b in zip(r_row, q_row)))
    return "".join(top), "".join(mid), "".join(bot)


# ---------------------------------------------------------------- analyzer
@dataclass
class MutationReport:
    reference: str
    query: str
    lcs: str
    lcs_positions: dict
    variants: list
    anchors: list
    gaps: list
    identity: float
    low_similarity: bool
    gst_nodes: int
    gst_text_length: int
    algorithm: str
    timings: dict = field(default_factory=dict)

    @property
    def counts(self) -> dict:
        out = {}
        for v in self.variants:
            out[v.type] = out.get(v.type, 0) + 1
        return out

    def alignment(self):
        return build_alignment(self.reference, self.query, self.anchors, self.gaps)


class MutationAnalyzer:
    def __init__(self, ref: str, qry: str, min_anchor: int = 3, algorithm: str = "naive",
                 flank: int = 4):
        if min_anchor < 1:
            raise ValueError("min_anchor must be >= 1")
        self.ref = validate_dna(ref, "reference")
        self.qry = validate_dna(qry, "query")
        self.min_anchor = min_anchor
        self.algorithm = algorithm
        self.flank = flank
        t0 = time.perf_counter()
        self.gst = GeneralizedSuffixTree([self.ref, self.qry], algorithm)
        self._build_time = time.perf_counter() - t0

    def report(self) -> MutationReport:
        timings = {"gst_construction": self._build_time}
        t0 = time.perf_counter()
        lcs, pos = self.gst.lcs_with_positions()
        timings["lcs"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        anchors, gaps = anchor_alignment(self.ref, self.qry, self.min_anchor, self.algorithm)
        variants = []
        for r0, r1, q0, q1 in gaps:
            variants.extend(classify_gap(self.ref[r0:r1], self.qry[q0:q1], r0 + 1))
        variants.sort(key=lambda v: (v.position, v.type))
        for v in variants:
            v.context = _context(self.ref, v, self.flank)
        timings["variant_detection"] = time.perf_counter() - t0

        matched = sum(a.length for a in anchors)
        identity = matched / max(len(self.ref), len(self.qry))
        return MutationReport(
            reference=self.ref, query=self.qry, lcs=lcs, lcs_positions=pos,
            variants=variants, anchors=anchors, gaps=gaps, identity=identity,
            low_similarity=identity < 0.5, gst_nodes=self.gst.count_nodes(),
            gst_text_length=len(self.gst.text), algorithm=self.algorithm, timings=timings)
