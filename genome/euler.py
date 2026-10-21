"""Part 3 - Eulerian traversal (Hierholzer), path spelling and unitigs.

Assembly = a walk that uses every k-mer (edge) exactly once = an Eulerian path.
Directed Euler conditions:
    cycle: in(v) == out(v) for every v
    path : exactly one start with out - in = +1, one end with in - out = +1,
           every other vertex balanced
plus all edges must lie in one connected component.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field


@dataclass
class EulerCheck:
    ok: bool
    kind: str                     # "path", "cycle", "empty" or "none"
    start: str = None
    end: str = None
    reason: str = ""
    unbalanced: dict = field(default_factory=dict)   # vertex -> out - in


def euler_endpoints(g):
    """(start, end, ok) from the degree conditions only. O(V)."""
    start = end = None
    for v in sorted(g.vertices):
        d = g.out_degree(v) - g.in_degree(v)
        if d == 1:
            if start is not None:
                return None, None, False
            start = v
        elif d == -1:
            if end is not None:
                return None, None, False
            end = v
        elif d != 0:
            return None, None, False
    if (start is None) != (end is None):
        return None, None, False
    return start, end, True


def check_euler(g) -> EulerCheck:
    """Degree conditions + connectivity, with a human-readable reason. O(V + E)."""
    if g.num_edges() == 0:
        return EulerCheck(False, "empty", reason="graph has no edges (no k-mers)")
    unbalanced = {v: g.out_degree(v) - g.in_degree(v) for v in sorted(g.vertices)
                  if g.out_degree(v) != g.in_degree(v)}
    start, end, ok = euler_endpoints(g)
    if not ok:
        starts = [v for v, d in unbalanced.items() if d > 0]
        ends = [v for v, d in unbalanced.items() if d < 0]
        return EulerCheck(False, "none", reason=(
            f"degree condition violated: {len(starts)} start candidate(s) {starts[:6]}, "
            f"{len(ends)} end candidate(s) {ends[:6]}"), unbalanced=unbalanced)
    comps = [c for c in g.weakly_connected_components()
             if any(g.out_degree(v) or g.in_degree(v) for v in c)]
    if len(comps) > 1:
        return EulerCheck(False, "none", reason=(
            f"degrees balance but the graph is split into {len(comps)} disconnected components"),
            unbalanced=unbalanced)
    if start is None:
        return EulerCheck(True, "cycle", reason="all vertices balanced: Eulerian cycle")
    return EulerCheck(True, "path", start, end,
                      reason=f"start {start} (out-in=+1), end {end} (in-out=+1)",
                      unbalanced=unbalanced)


def hierholzer(g, start: str = None):
    """Iterative Hierholzer. Returns the vertex list of an Euler path/cycle, or None.

    Each vertex keeps a pointer to its next unused edge, so an edge is consumed in
    O(1) without deleting it. Every edge is pushed and popped once: O(E) time,
    O(V + E) space.
    """
    s, _end, ok = euler_endpoints(g)
    if not ok:
        return None
    if start is None:
        start = s
    if start is None:                                   # cycle: any vertex with edges
        start = next((v for v in sorted(g.vertices) if g.out_degree(v) > 0), None)
        if start is None:
            return []
    ptr = {v: 0 for v in g.vertices}
    stack, path = [start], []
    while stack:
        v = stack[-1]
        if ptr[v] < len(g.out[v]):
            e = g.out[v][ptr[v]]
            ptr[v] += 1                                 # consume the edge
            stack.append(e.dst)
        else:
            path.append(stack.pop())                    # stuck: v is finished
    path.reverse()
    if len(path) != g.num_edges() + 1:                  # some edges unreachable
        return None
    return path


def path_to_sequence(path) -> str:
    """First vertex in full, then the last character of every following vertex."""
    if not path:
        return ""
    return path[0] + "".join(v[-1] for v in path[1:])


def unitigs(g) -> list:
    """Maximal non-branching paths (safe contigs). Isolated cycles are included. O(V + E)."""
    def one_in_one_out(v):
        return g.in_degree(v) == 1 and g.out_degree(v) == 1

    used, result = set(), []
    for v in sorted(g.vertices):
        if one_in_one_out(v):
            continue
        for e in g.out[v]:
            used.add(id(e))
            seq, cur = v + e.dst[-1], e.dst
            while one_in_one_out(cur):
                nxt = g.out[cur][0]
                if id(nxt) in used:
                    break
                used.add(id(nxt))
                seq += nxt.dst[-1]
                cur = nxt.dst
            result.append(seq)
    for v in sorted(g.vertices):                        # leftover isolated cycles
        for e in g.out[v]:
            if id(e) in used:
                continue
            seq, cur_e = v, e
            while cur_e is not None and id(cur_e) not in used:
                used.add(id(cur_e))
                seq += cur_e.dst[-1]
                nxt = g.out[cur_e.dst]
                cur_e = nxt[0] if len(nxt) == 1 else None
            result.append(seq)
    return result


def enumerate_euler_sequences(g, limit: int = 1000, max_steps: int = 200_000):
    """Brute-force distinct sequences spelled by Euler paths (tiny graphs only).

    Returns (sequences, exhaustive). ``exhaustive`` is False if the search was cut
    off by ``max_steps`` before the space was fully explored.
    """
    start, _end, ok = euler_endpoints(g)
    if not ok or g.num_edges() == 0:
        return set(), True
    starts = [start] if start else sorted(v for v in g.vertices if g.out_degree(v))
    m = g.num_edges()
    used = {id(e): False for e in g.edges()}
    found = set()
    steps = [0]
    exhaustive = [True]
    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old_limit, m + 200))

    def dfs(v, seq, n):
        steps[0] += 1
        if steps[0] > max_steps:
            exhaustive[0] = False
            return
        if len(found) >= limit:
            return
        if n == m:
            found.add(seq)
            return
        tried = set()
        for e in g.out[v]:
            if not used[id(e)] and e.dst not in tried:
                tried.add(e.dst)
                used[id(e)] = True
                dfs(e.dst, seq + e.dst[-1], n + 1)
                used[id(e)] = False

    try:
        for s0 in starts:
            dfs(s0, s0, 0)
    finally:
        sys.setrecursionlimit(old_limit)
    if len(found) >= limit:
        exhaustive[0] = False
    return found, exhaustive[0]


def all_euler_sequences(g, limit: int = 1000) -> set:
    return enumerate_euler_sequences(g, limit)[0]
