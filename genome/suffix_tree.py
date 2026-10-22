"""Part 4 - Generalized Suffix Tree (GST).

Several strings are indexed in ONE compressed suffix tree. Each string gets its
own unique sentinel, the strings are concatenated (ref + '#' + query + '$') and
every suffix of that text is inserted. Edge labels are stored as (start, end)
index pairs into the text, so the tree has O(n) nodes.

Two construction algorithms are provided:
    "naive"   - insert every suffix from the root, splitting edges where it
                diverges. O(n^2) worst case (e.g. AAAA...), fast on random DNA.
    "ukkonen" - Ukkonen's online algorithm with suffix links. O(n) for a
                constant alphabet.
Both produce exactly the same tree.
"""

from __future__ import annotations

SENTINELS = "#$%&!@^*"


class Node:
    __slots__ = ("children", "start", "end", "link", "leaf_info", "depth")

    def __init__(self, start: int = -1, end=-1):
        self.children = {}      # first character of edge -> Node
        self.start = start      # edge label into this node = text[start:end]
        self.end = end          # None while building (Ukkonen leaf: "global end")
        self.link = None        # suffix link (Ukkonen only)
        self.leaf_info = None   # (sequence id, offset) for leaves
        self.depth = 0          # length of the path label root -> this node

    @property
    def is_leaf(self) -> bool:
        return not self.children


class GeneralizedSuffixTree:
    def __init__(self, seqs, algorithm: str = "naive"):
        if isinstance(seqs, str):
            seqs = [seqs]
        seqs = list(seqs)
        if not seqs:
            raise ValueError("need at least one sequence")
        if len(seqs) > len(SENTINELS):
            raise ValueError(f"at most {len(SENTINELS)} sequences are supported")
        for s in seqs:
            bad = set(s) & set(SENTINELS)
            if bad:
                raise ValueError(f"sequence contains reserved sentinel character(s) {sorted(bad)}")
        if algorithm not in ("naive", "ukkonen"):
            raise ValueError("algorithm must be 'naive' or 'ukkonen'")
        self.seqs = seqs
        self.algorithm = algorithm
        self.sentinels = SENTINELS[:len(seqs)]
        parts, owner = [], []
        for i, s in enumerate(seqs):
            parts.append(s + self.sentinels[i])
            owner.extend((i, j) for j in range(len(s) + 1))
        self.text = "".join(parts)
        self.owner = owner              # text position -> (sequence id, offset)
        self.root = Node()
        self._gend = 0
        if algorithm == "naive":
            self._build_naive()
        else:
            self._build_ukkonen()
        self._finalize()
        self._stats = None

    # ------------------------------------------------------------ construction
    def _build_naive(self):
        t, n = self.text, len(self.text)
        for i in range(n):                              # insert suffix t[i:]
            node, pos = self.root, i
            while True:
                child = node.children.get(t[pos])
                if child is None:
                    node.children[t[pos]] = Node(pos, n)
                    break
                j = child.start                         # walk down the edge
                while j < child.end and pos < n and t[j] == t[pos]:
                    j += 1
                    pos += 1
                if j == child.end:
                    node = child
                    continue
                mid = Node(child.start, j)              # split the edge at the mismatch
                node.children[t[child.start]] = mid
                child.start = j
                mid.children[t[j]] = child
                mid.children[t[pos]] = Node(pos, n)
                break

    def _edge_len(self, node: Node) -> int:
        return (self._gend if node.end is None else node.end) - node.start

    def _build_ukkonen(self):
        t, n, root = self.text, len(self.text), self.root
        active_node, active_edge, active_length, remaining = root, 0, 0, 0
        for i in range(n):
            self._gend = i + 1                          # rule 1: every leaf grows by one
            remaining += 1
            last_new = None
            c = t[i]
            while remaining > 0:
                if active_length == 0:
                    active_edge = i
                nxt = active_node.children.get(t[active_edge])
                if nxt is None:                         # rule 2: new leaf at active node
                    active_node.children[t[active_edge]] = Node(i, None)
                    if last_new is not None:
                        last_new.link = active_node
                        last_new = None
                else:
                    el = self._edge_len(nxt)
                    if active_length >= el:             # skip/count trick
                        active_edge += el
                        active_length -= el
                        active_node = nxt
                        continue
                    if t[nxt.start + active_length] == c:   # rule 3: already present
                        if last_new is not None and active_node is not root:
                            last_new.link = active_node
                            last_new = None
                        active_length += 1
                        break
                    split = Node(nxt.start, nxt.start + active_length)   # rule 2: split
                    split.link = root
                    active_node.children[t[active_edge]] = split
                    split.children[c] = Node(i, None)
                    nxt.start += active_length
                    split.children[t[nxt.start]] = nxt
                    if last_new is not None:
                        last_new.link = split
                    last_new = split
                remaining -= 1
                if active_node is root and active_length > 0:
                    active_length -= 1
                    active_edge = i - remaining + 1
                elif active_node is not root:
                    active_node = active_node.link if active_node.link is not None else root

    def _finalize(self):
        """Close open leaf ends, compute depths and label every leaf (seq id, offset)."""
        n = len(self.text)
        stack = [(self.root, 0)]
        while stack:
            node, depth = stack.pop()
            node.depth = depth
            if node.end is None:
                node.end = n
            if not node.children and node is not self.root:
                node.leaf_info = self.owner[n - depth]
            for child in node.children.values():
                if child.end is None:
                    child.end = n
                stack.append((child, depth + child.end - child.start))

    # ------------------------------------------------------------ queries
    def _locate(self, pattern: str):
        """Node at or below the end of ``pattern``'s path, or None. O(m)."""
        node, i, t = self.root, 0, self.text
        while i < len(pattern):
            child = node.children.get(pattern[i])
            if child is None:
                return None
            j = child.start
            while j < child.end and i < len(pattern):
                if t[j] != pattern[i]:
                    return None
                j += 1
                i += 1
            node = child
        return node

    def contains(self, pattern: str) -> bool:
        if not pattern:
            return True
        return self._locate(pattern) is not None

    def occurrences(self, pattern: str) -> list:
        """All (sequence id, offset) where ``pattern`` occurs. O(m + occ)."""
        if not pattern or set(pattern) & set(self.sentinels):
            return []
        node = self._locate(pattern)
        if node is None:
            return []
        out, stack = [], [node]
        while stack:
            cur = stack.pop()
            if cur.leaf_info is not None:
                out.append(cur.leaf_info)
            stack.extend(cur.children.values())
        return sorted(out)

    def count(self, pattern: str) -> int:
        return len(self.occurrences(pattern))

    def _node_stats(self):
        """Post-order pass: per node, leaf count and first offset per sequence. O(n)."""
        if self._stats is not None:
            return self._stats
        m = len(self.seqs)
        order, stack = [], [self.root]
        while stack:
            node = stack.pop()
            order.append(node)
            stack.extend(node.children.values())
        stats = {}
        for node in reversed(order):                    # children before parents
            counts, first = [0] * m, [None] * m
            if node.leaf_info is not None:
                sid, off = node.leaf_info
                counts[sid], first[sid] = 1, off
            for child in node.children.values():
                c_counts, c_first = stats[id(child)]
                for s in range(m):
                    counts[s] += c_counts[s]
                    if c_first[s] is not None and (first[s] is None or c_first[s] < first[s]):
                        first[s] = c_first[s]
            stats[id(node)] = (counts, first)
        self._stats = (order, stats)
        return self._stats

    def _path_label(self, node: Node) -> str:
        """Internal nodes never contain a sentinel, so this is a plain substring."""
        _, stats = self._node_stats()
        counts, first = stats[id(node)]
        sid = next(s for s in range(len(self.seqs)) if first[s] is not None)
        return self.seqs[sid][first[sid]:first[sid] + node.depth]

    def lcs_with_positions(self, seq_ids=None):
        """Longest common substring of the given sequences (default: all).

        The deepest internal node whose subtree has leaves from every sequence.
        Ties are broken by the earliest offset in the first sequence, then the
        second. Returns (substring, {seq_id: first offset}) or ("", {}).
        """
        seq_ids = list(range(len(self.seqs))) if seq_ids is None else list(seq_ids)
        order, stats = self._node_stats()
        best_key, best_node = None, None
        for node in order:
            if node is self.root or node.leaf_info is not None or node.depth == 0:
                continue
            counts, first = stats[id(node)]
            if all(counts[s] > 0 for s in seq_ids):
                key = (node.depth, tuple(-first[s] for s in seq_ids))
                if best_key is None or key > best_key:
                    best_key, best_node = key, node
        if best_node is None:
            return "", {}
        _, first = stats[id(best_node)]
        return self._path_label(best_node), {s: first[s] for s in seq_ids}

    def longest_common_substring(self, seq_ids=None) -> str:
        return self.lcs_with_positions(seq_ids)[0]

    def longest_repeated_substring(self, seq_id: int = 0) -> str:
        """Deepest internal node with >= 2 leaves from the same sequence."""
        order, stats = self._node_stats()
        best = None
        for node in order:
            if node is self.root or node.leaf_info is not None:
                continue
            counts, _ = stats[id(node)]
            if counts[seq_id] >= 2 and (best is None or node.depth > best.depth):
                best = node
        if best is None:
            return ""
        _, first = stats[id(best)]
        return self.seqs[seq_id][first[seq_id]:first[seq_id] + best.depth]

    def count_nodes(self) -> int:
        return len(self._node_stats()[0])

    def count_leaves(self) -> int:
        return sum(1 for n in self._node_stats()[0] if n.leaf_info is not None)

    def edge_label(self, node: Node) -> str:
        """Edge label, cut after the first sentinel (a real GST suffix ends at ITS sentinel)."""
        label = self.text[node.start:node.end]
        cuts = [label.index(s) for s in self.sentinels if s in label]
        return label[:min(cuts) + 1] if cuts else label

    def dump(self, max_lines: int = 400) -> str:
        """Indented text drawing of the tree (children in sorted order)."""
        lines = ["root"]
        stack = [(c, 1) for _, c in sorted(self.root.children.items(), reverse=True)]
        while stack:
            if len(lines) >= max_lines:
                lines.append(f"... (truncated, {self.count_nodes()} nodes in total)")
                break
            node, indent = stack.pop()
            tag = ""
            if node.leaf_info is not None:
                sid, off = node.leaf_info
                tag = f"   <- leaf (seq{sid + 1}, pos {off})"
            lines.append("  " * indent + "+-- " + self.edge_label(node) + tag)
            stack.extend((c, indent + 1) for _, c in sorted(node.children.items(), reverse=True))
        return "\n".join(lines)

    def signature(self) -> list:
        """Sorted (path label, leaf info) pairs - used to prove naive == Ukkonen."""
        out, stack = [], [(self.root, "")]
        while stack:
            node, label = stack.pop()
            out.append((label, node.leaf_info))
            for child in node.children.values():
                stack.append((child, label + self.text[child.start:child.end]))
        return sorted(out, key=lambda x: (x[0], x[1] or (-1, -1)))


class SuffixTree(GeneralizedSuffixTree):
    """Suffix tree of a single string (string + '#')."""

    def __init__(self, text: str, algorithm: str = "naive"):
        super().__init__([text], algorithm)
