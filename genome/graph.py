"""Part 2 - Graph ADT implemented as a De Bruijn graph (adjacency lists).

Convention (edge-centric, as in Pevzner/Compeau): vertices are (k-1)-mers and
every distinct k-mer is a directed edge prefix -> suffix carrying its coverage.

Representation: two hash maps of dynamic arrays
    out[v] = list of edges leaving v      (adjacency list)
    inc[v] = list of edges entering v     (reverse adjacency list)
A vertex has at most 4 outgoing edges (next base A/C/G/T), so E <= 4V and the
graph is extremely sparse; an adjacency matrix would need (4^(k-1))^2 cells.
"""

from __future__ import annotations

from collections import deque


class Edge:
    """One k-mer: src (prefix) -> dst (suffix) with its observed coverage."""

    __slots__ = ("src", "dst", "kmer", "coverage")

    def __init__(self, src: str, dst: str, kmer: str, coverage: int = 1):
        self.src, self.dst, self.kmer, self.coverage = src, dst, kmer, coverage

    def __repr__(self):
        return f"{self.src}->{self.dst}[{self.kmer} x{self.coverage}]"


def _mean_coverage(edges) -> float:
    return sum(e.coverage for e in edges) / len(edges)


class DeBruijnGraph:
    def __init__(self, k: int):
        if k < 2:
            raise ValueError("k must be >= 2 (vertices are (k-1)-mers)")
        self.k = k
        self.out = {}   # vertex -> [Edge]
        self.inc = {}   # vertex -> [Edge]
        self._num_edges = 0

    # ------------------------------------------------------------ Graph ADT
    @property
    def vertices(self):
        return self.out.keys()

    def add_vertex(self, v: str) -> None:
        if v not in self.out:
            self.out[v] = []
            self.inc[v] = []

    def add_edge(self, src: str, dst: str, kmer: str = None, coverage: int = 1) -> Edge:
        """Insert edge src -> dst. O(k) for the overlap check, O(1) amortised insert."""
        if len(src) != self.k - 1 or len(dst) != self.k - 1:
            raise ValueError(f"vertices must be (k-1)-mers of length {self.k - 1}")
        if src[1:] != dst[:-1]:
            raise ValueError(f"{src} and {dst} do not overlap by k-2 characters")
        if kmer is None:
            kmer = src + dst[-1]
        self.add_vertex(src)
        self.add_vertex(dst)
        edge = Edge(src, dst, kmer, coverage)
        self.out[src].append(edge)
        self.inc[dst].append(edge)
        self._num_edges += 1
        return edge

    def remove_edge(self, edge: Edge) -> None:
        """O(deg) = O(4): list.remove on a tiny adjacency list."""
        self.out[edge.src].remove(edge)
        self.inc[edge.dst].remove(edge)
        self._num_edges -= 1

    def remove_vertex(self, v: str) -> None:
        for e in list(self.out[v]) + list(self.inc[v]):
            if e in self.out.get(e.src, ()):
                self.remove_edge(e)
        del self.out[v]
        del self.inc[v]

    def remove_isolated_vertices(self) -> int:
        isolated = [v for v in self.out if not self.out[v] and not self.inc[v]]
        for v in isolated:
            del self.out[v]
            del self.inc[v]
        return len(isolated)

    def neighbors(self, v: str) -> list:
        return [e.dst for e in self.out[v]]

    def out_degree(self, v: str) -> int:
        return len(self.out[v])

    def in_degree(self, v: str) -> int:
        return len(self.inc[v])

    def has_vertex(self, v: str) -> bool:
        return v in self.out

    def has_edge(self, u: str, v: str) -> bool:
        """O(outdeg(u)) <= O(4)."""
        return u in self.out and any(e.dst == v for e in self.out[u])

    def edges(self):
        for v in self.out:
            yield from self.out[v]

    def num_vertices(self) -> int:
        return len(self.out)

    def num_edges(self) -> int:
        return self._num_edges

    def is_empty(self) -> bool:
        return not self.out

    # ------------------------------------------------------------ construction
    @classmethod
    def from_kmer_counts(cls, freq, k: int, min_coverage: int = 1, multigraph: bool = False):
        """Build the graph from a {k-mer: count} table.

        multigraph=False: one edge per DISTINCT k-mer, count kept as coverage (reads).
        multigraph=True : one edge per OCCURRENCE (freq = multiplicity inside one
                          known sequence, so genuine repeats are kept).
        k-mers seen fewer than ``min_coverage`` times are dropped (error filter).
        Time O(D * k) for D distinct k-mers; space O(V + E).
        """
        g = cls(k)
        for kmer, count in sorted(freq.items()):
            if len(kmer) != k:
                raise ValueError(f"k-mer {kmer!r} does not have length {k}")
            if count >= min_coverage:
                for _ in range(count if multigraph else 1):
                    g.add_edge(kmer[:-1], kmer[1:], kmer, count)
        return g

    @classmethod
    def from_sequence(cls, seq: str, k: int, multigraph: bool = True):
        from collections import Counter
        freq = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
        return cls.from_kmer_counts(freq, k, multigraph=multigraph)

    def copy(self):
        g = DeBruijnGraph(self.k)
        for v in self.out:
            g.add_vertex(v)
        for e in self.edges():
            g.add_edge(e.src, e.dst, e.kmer, e.coverage)
        return g

    # ------------------------------------------------------------ queries
    def degree_report(self) -> dict:
        return {v: (self.in_degree(v), self.out_degree(v)) for v in sorted(self.out)}

    def branching_vertices(self) -> list:
        return sorted(v for v in self.out if self.out_degree(v) > 1 or self.in_degree(v) > 1)

    def weakly_connected_components(self) -> list:
        """BFS ignoring edge direction. O(V + E)."""
        seen, comps = set(), []
        for s in sorted(self.out):
            if s in seen:
                continue
            comp, queue = set(), deque([s])
            seen.add(s)
            while queue:
                v = queue.popleft()
                comp.add(v)
                for nb in [e.dst for e in self.out[v]] + [e.src for e in self.inc[v]]:
                    if nb not in seen:
                        seen.add(nb)
                        queue.append(nb)
            comps.append(comp)
        return comps

    def summary(self) -> dict:
        return {
            "k": self.k,
            "vertices": self.num_vertices(),
            "edges": self.num_edges(),
            "branching_vertices": len(self.branching_vertices()),
            "components": len(self.weakly_connected_components()),
        }

    # ------------------------------------------------------------ cleaning
    def filter_low_coverage(self, min_coverage: int) -> int:
        """Remove edges observed fewer than ``min_coverage`` times. O(E)."""
        doomed = [e for e in self.edges() if e.coverage < min_coverage]
        for e in doomed:
            self.remove_edge(e)
        self.remove_isolated_vertices()
        return len(doomed)

    def _tip_chain(self, v: str, max_length: int):
        """If v is the loose end of a tip, return (chain_edges, alternative_coverage)."""
        if self.out_degree(v) == 0 and self.in_degree(v) == 1:
            # dead end: walk backwards to the branching vertex it hangs from
            edges, cur = [], v
            while True:
                e = self.inc[cur][0]
                edges.append(e)
                cur = e.src
                if self.out_degree(cur) >= 2:
                    alt = max(x.coverage for x in self.out[cur] if x is not e)
                    return edges, alt
                if (len(edges) >= max_length or cur == v
                        or self.in_degree(cur) != 1 or self.out_degree(cur) != 1):
                    return None
        if self.in_degree(v) == 0 and self.out_degree(v) == 1:
            # dangling start: walk forwards to the vertex where it joins
            edges, cur = [], v
            while True:
                e = self.out[cur][0]
                edges.append(e)
                cur = e.dst
                if self.in_degree(cur) >= 2:
                    alt = max(x.coverage for x in self.inc[cur] if x is not e)
                    return edges, alt
                if (len(edges) >= max_length or cur == v
                        or self.in_degree(cur) != 1 or self.out_degree(cur) != 1):
                    return None
        return None

    def remove_tips(self, max_length: int = None, max_ratio: float = 0.5) -> int:
        """Delete short dead-end branches (typically errors near read ends).

        A tip is removed only when it is at most ``max_length`` edges long
        (default 2k) and its mean coverage is <= max_ratio x the coverage of the
        competing branch, so genuine sequence ends with normal coverage survive.
        """
        max_length = max_length or 2 * self.k
        removed, changed = 0, True
        while changed:
            changed = False
            for v in sorted(self.out):
                if v not in self.out:
                    continue
                tip = self._tip_chain(v, max_length)
                if tip is None:
                    continue
                edges, alt = tip
                if _mean_coverage(edges) <= max_ratio * alt:
                    for e in edges:
                        self.remove_edge(e)
                    removed += len(edges)
                    changed = True
            self.remove_isolated_vertices()
        return removed

    def _walk_chain(self, first: Edge, max_length: int):
        """Follow non-branching (in=1, out=1) vertices starting with edge ``first``."""
        edges, cur = [first], first.dst
        while (self.in_degree(cur) == 1 and self.out_degree(cur) == 1
               and len(edges) < max_length and cur != first.src):
            e = self.out[cur][0]
            edges.append(e)
            cur = e.dst
        return edges, cur

    def _distance(self, start: str, target: str, max_depth: int):
        """BFS distance (in edges) from start to target, or None if > max_depth."""
        if start == target:
            return 0
        seen, frontier = {start}, [start]
        for depth in range(1, max_depth + 1):
            nxt = []
            for v in frontier:
                for e in self.out[v]:
                    if e.dst == target:
                        return depth
                    if e.dst not in seen:
                        seen.add(e.dst)
                        nxt.append(e.dst)
            frontier = nxt
        return None

    def pop_bubbles(self, max_length: int = None, max_ratio: float = 0.5) -> int:
        """Collapse bubbles: two paths that leave the same vertex and meet again.

        For each branching vertex u, each outgoing non-branching chain is walked to
        the vertex w where it merges (in-degree >= 2). If another branch of u also
        reaches w by a path of about the same length (+-2 edges: a substitution
        error gives two branches of exactly k edges) and has at least 1/max_ratio
        times the coverage, the weaker chain is an error path and is deleted.
        """
        max_length = max_length or 2 * self.k
        removed, changed = 0, True
        while changed:
            changed = False
            for u in sorted(self.out):
                if u not in self.out or self.out_degree(u) < 2:
                    continue
                branches = [self._walk_chain(e, max_length) for e in self.out[u]]
                branches.sort(key=lambda b: _mean_coverage(b[0]))
                for i, (edges, end) in enumerate(branches):
                    if end == u or self.in_degree(end) < 2:
                        continue
                    weak = _mean_coverage(edges)
                    for j, (other, _) in enumerate(branches):
                        if j == i or weak > max_ratio * _mean_coverage(other):
                            continue
                        d = self._distance(other[0].dst, end, len(edges) + 2)
                        if d is not None and abs(1 + d - len(edges)) <= 2:
                            for e in edges:
                                self.remove_edge(e)
                            removed += len(edges)
                            changed = True
                            break
                    if changed:
                        break
            self.remove_isolated_vertices()
        return removed

    # ------------------------------------------------------------ export
    def to_dot(self, path: list = None) -> str:
        """Graphviz DOT text; with an Euler ``path`` the edges are numbered in walk order."""
        order = {}
        if path:
            for step, (a, b) in enumerate(zip(path, path[1:]), start=1):
                order.setdefault((a, b), []).append(str(step))
        branch = set(self.branching_vertices())
        lines = ["digraph DeBruijn {", "  rankdir=LR;", "  node [shape=circle, fontname=Helvetica];"]
        for v in sorted(self.out):
            style = ' style=filled fillcolor="#fff3c4"' if v in branch else ""
            lines.append(f'  "{v}" [label="{v}"{style}];')
        for e in sorted(self.edges(), key=lambda e: (e.src, e.kmer)):
            label = f"{e.kmer} x{e.coverage}"
            steps = order.get((e.src, e.dst))
            if steps:
                label += f" [#{steps.pop(0)}]"
            lines.append(f'  "{e.src}" -> "{e.dst}" [label="{label}"];')
        lines.append("}")
        return "\n".join(lines) + "\n"
