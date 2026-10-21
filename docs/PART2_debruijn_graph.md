# Part 2 — Graph ADT: the De Bruijn Graph

**Files:** `genome/graph.py`, `gui.py`
**Tests:** `tests/test_part2_graph.py`
**Guide sections:** Parts 4, 5, 9, 11

## What this part does
It builds the assembly graph from the k-mer table. **Vertices are (k−1)-mers and every distinct
k-mer is a directed edge prefix → suffix that carries its coverage.** It provides the Graph ADT
operations (add/remove vertex and edge, degrees, neighbours, components) and the graph-cleaning
operations that deal with sequencing errors: coverage filter, tip removal and bubble popping.
It also draws the graph in the GUI.

## 1. Data structure justification
**Why a graph?** Assembly asks "which piece can follow which?". That is a pairwise *relationship*
between overlapping pieces, which is exactly what edges model.

**Why k-mers as edges (not vertices)?** "Use every k-mer exactly once" then becomes "traverse
every edge exactly once", an **Eulerian path**, which is solvable in linear time. If reads or
k-mers were vertices we would need a **Hamiltonian path**, which is NP-complete.

**Why (k−1)-mers as vertices?** Two consecutive k-mers (ATG, TGC) overlap in exactly k−1
letters (TG). That shared piece is where the two edges meet, so it has to be the vertex.

**Adjacency list vs adjacency matrix:**

| | Adjacency list (chosen) | Adjacency matrix |
|---|---|---|
| Space | O(V + E) | O(V²) |
| Neighbours of v | O(outdeg v) ≤ O(4) | O(V) |
| hasEdge(u,v) | O(outdeg u) ≤ O(4) | O(1) |
| Add vertex | O(1) amortised | O(V²) resize |

Every vertex has at most 4 outgoing edges (next base A/C/G/T), so **E ≤ 4V**: the graph is
extremely sparse. The number of *possible* vertices is 4^(k−1); for k = 31 that is about
1.15×10¹⁸, so a matrix is physically impossible. Even V = 10⁶ would need 10¹² cells. The list
stores only the vertices that actually occur.

## 2. Implementation strategy (code-level choices)
- `out: dict[vertex → list[Edge]]` and `inc: dict[vertex → list[Edge]]`. Each is a **hash map of dynamic arrays**. The hash map finds a vertex by its string label in O(1) expected time; the small array holds at most 4 edges.
- We keep a **reverse adjacency list (`inc`)** in addition to `out`. It costs 2× edge references, but in-degree becomes O(1) and tip removal can walk *backwards* from a dead end. Without it, we would have to scan every edge to find a vertex's predecessors (O(E)).
- `Edge` uses `__slots__`, so each edge has no per-object `__dict__`, which matters with thousands of edges.
- **One edge per distinct k-mer**, with its count stored as `coverage`. Overlapping reads therefore don't create parallel edges that would break the degree balance. `multigraph=True` (one edge per occurrence) is available for a single known sequence, where repeats must be kept (guide Part 34).
- `add_edge` checks that src and dst really overlap by k−2 characters, so an invalid graph can't be built.
- `num_edges` is a maintained counter, O(1), used by Hierholzer's "all edges used?" check.

**Cleaning (handling errors):**
- `filter_low_coverage(t)`: delete k-mers seen fewer than t times. A sequencing error creates up to k fake k-mers, each with coverage ≈ 1.
- `remove_tips()`: a tip is a short dead end (≤ 2k edges) hanging off a branching vertex. It is deleted only if its coverage is ≤ 50% of the competing branch, so a real sequence end with normal coverage survives.
- `pop_bubbles()`: two branches leave vertex u and meet again at w with the **same length (±2 edges)**. If one has ≤ 50% of the other's coverage, it is an error path and is removed. A *balanced* bubble (two real alleles) is kept. The same-length rule exists because, without it, a genuine edge next to a high-coverage detour was being deleted (see `test_balanced_bubble_is_kept`).

## 3. Time complexity (with justification)
- **Construction O(D·k)**, D = distinct k-mers: each k-mer does two O(k) slices (prefix, suffix), two O(1)-expected hash lookups and two O(1) amortised appends.
- **Space O(V + E)**, E = D. For a connected graph V ≤ E + 1.
- **add_edge O(k)** (the overlap check), **remove_edge O(deg) = O(4)** (`list.remove` on a tiny list), **degrees O(1)**, **has_edge O(4)**.
- **Coverage filter O(E)**: one pass over the edges.
- **Components O(V + E)**: BFS visits every vertex and edge once.
- **Tip removal O(V·k)**: every vertex is a candidate and each walk stops after ≤ 2k edges.
- **Bubble popping**: per branching vertex, each branch walk is ≤ 2k edges and the meeting check is a BFS limited to the branch length + 2. This bounded local search is linear in V for fixed k.

## 4. Live demo (my part)
```bash
python main.py assemble --reads data/reads.fasta -k 3                 # graph listing: AT has out-degree 2 (branch)
python main.py assemble --reads data/error_reads.fasta -k 3           # error -> bubble CGT/GTT/TTA (coverage 1)
python main.py assemble --reads data/error_reads.fasta -k 3 --clean   # bubble popped -> ATGCGATAC
python main.py assemble --reads data/reads.fasta -k 3 --out out       # out/graph.dot for Graphviz
python main.py gui                                                    # drawn graph, Euler order numbered
```

## 5. Questions I should be ready for
- *Why not an adjacency matrix?* E ≤ 4V, and 4^(k−1) possible vertices make it O(V²) and impossible.
- *What is a branch / tip / bubble?* out-degree ≥ 2 / a short dead end (error near a read end) / two paths between the same vertices (error, or a real variant if both have similar coverage).
- *Why store coverage on the edge?* It is the evidence used to tell errors from real sequence.
- *What happens to duplicate k-mers?* With the distinct-edge model they collapse into one edge and the count goes up. Genuine repeats inside one genome are therefore collapsed too, unless you use `multigraph`.
- *Why `inc` as well as `out`?* O(1) in-degree and backward walks, at the cost of 2× edge references.

## Suggested commit order
1. `Edge` + `DeBruijnGraph` core ADT (add/remove/degree/neighbours)
2. `from_kmer_counts` / `from_sequence` construction
3. Components, summary, degree report, DOT export
4. Cleaning: coverage filter, tips, bubbles
5. `tests/test_part2_graph.py`
6. `gui.py` (after the pipeline is merged)
