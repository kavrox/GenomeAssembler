# Part 3 — Eulerian Traversal & the Assembly Pipeline

**Files:** `genome/euler.py`, `genome/assembler.py`, `genome/simulate.py`, `genome/benchmark.py`, `demo.py`
**Tests:** `tests/test_part3_euler_assembly.py`
**Guide sections:** Parts 6, 7, 8, 10, 21

## What this part does
It checks Euler's conditions on the De Bruijn graph and finds an Eulerian path with **Hierholzer's
algorithm**. It then spells the path into the assembled sequence. If no Euler path exists, it
falls back to **unitigs** (maximal non-branching paths). It also runs the whole assembly pipeline,
detects repeat ambiguity, runs the k-sweep experiment, simulates realistic data and measures Big-O
empirically.

## 1. Data structure / algorithm justification
**Eulerian, not Hamiltonian.** k-mers are edges, and a correct assembly must use every observed
k-mer exactly once. That means traversing every edge once: an Euler path. Testing and finding one
is **linear** (Hierholzer, O(E)). The alternative formulation (reads as vertices, visit every
vertex once) is a Hamiltonian path, which is **NP-complete**.

**Euler conditions (directed):**
- cycle: in(v) = out(v) for every v
- path: exactly one start with out − in = +1, one end with in − out = +1, all others balanced
- plus connectivity: all edges must be in one component. Balanced degrees alone are not enough; `test_balanced_but_disconnected` shows two separate cycles failing.

**Stack for Hierholzer.** "Go back to the latest vertex that still has unused edges" is LIFO,
which is exactly a stack (push O(1), pop O(1)). A recursive version would hit Python's recursion
limit around 1000 edges; the iterative stack works for any size.

**Edge pointer instead of deleting edges.** `ptr[v]` is the index of v's next unused edge.
Consuming an edge is `ptr[v] += 1`, which is O(1) and leaves the graph untouched. Physically
removing edges would be O(deg) per removal and would destroy the graph for later use (drawing,
reports).

## 2. Implementation strategy
- `check_euler` returns a reason ("degree condition violated: 2 start candidates [...]" / "split into 2 components"), so the demo explains *why* assembly failed.
- `hierholzer`: start at the vertex with out − in = +1 (or, for a cycle, the first vertex in sorted order). Walk unused edges, pushing as we go; when a vertex has no unused edge, pop it into `path`; reverse at the end. **Check `len(path) == E + 1`**: if not, some edges were never reachable and the graph is disconnected.
- `path_to_sequence`: write the first vertex in full, then add only the **last letter** of each following vertex, since consecutive vertices overlap by k−2. Length = (k−1) + E.
- Ties (several unused edges) are taken in adjacency order, which is sorted by k-mer, so output is deterministic.
- `unitigs`: start at every vertex that is *not* 1-in-1-out and extend through 1-in-1-out vertices. Pure isolated cycles are added at the end.
- **Ambiguity check**: for small graphs (≤ 400 edges), a bounded brute-force DFS looks for a *second* distinct Euler sequence. If one exists, the result is flagged as ambiguous (a repeat).
- `assemble()` times every stage (`perf_counter`) for the complexity slide.

## 3. Time complexity (with justification)
- **Degree check O(V)**: one pass computing out − in per vertex. **Connectivity O(V + E)**: BFS.
- **Hierholzer O(E)**: `ptr[v]` only moves forward, so each edge is consumed exactly once. Each consumption pushes one vertex, and each vertex is popped exactly once per push. Total pushes = E + 1, pops = E + 1, so the work is O(E) plus O(V) to initialise `ptr`.
- **Spelling O(V + E)**, i.e. the length of the path.
- **Unitigs O(V + E)**: each edge is added to exactly one unitig (tracked by the `used` set).
- **Whole assembly pipeline:** O(K·k) counting + O(D·k) graph + O(V+E) Euler = **linear in the input size**. `python main.py benchmark` confirms it: each stage takes about 2× as long when the genome doubles.
- Ambiguity enumeration is exponential in the worst case, which is why it is capped (≤ 400 edges, 50 000 steps). It is only a diagnostic.

## 4. Live demo (my part)
```bash
python demo.py        # k=3: AT -> TG -> GC -> CG -> GA -> AT -> TA -> AC  => ATGCGATAC
python demo.py 4      # k too large: degree condition fails, unitigs ATGC + GCGATAC
python demo.py 2      # k too small: valid Euler path spells TACGATGC (wrong) + ambiguity flagged
python main.py ksweep --reads data/reads.fasta --k 2-5
python main.py assemble --reads data/repeat.fasta -k 3 --multigraph     # TTGTGAGTA: 2 Euler paths
python main.py benchmark --quick
```

## 5. Questions I should be ready for
- *Why Eulerian and not Hamiltonian?* Edge-once is linear; vertex-once is NP-complete.
- *Trace Hierholzer on the Part 18 graph.* It first takes ATA → TAC, gets stuck at AC and pops it, then takes ATG … GAT back to AT. Reversing gives the path. (Guide Part 8.2; the code produces exactly this path.)
- *Why reverse the path?* Vertices are popped when finished, so dead ends (the true end) come out first.
- *What if two Euler paths exist?* That is repeat ambiguity (TTGTGAGTA vs TTGAGTGTA at k=3). Both are consistent with the data. A larger k resolves it (k=4 gives a unique path).
- *What if no Euler path exists?* Errors, low coverage or k too large. We report unitigs, the "safe" pieces every assembly must contain.
- *How does k affect the result?* Small k gives false overlaps (k=2 assembles the wrong string). Large k gives too little overlap, and the graph fragments (k=4). Choose the largest k your read overlap supports.
- *Path vs cycle?* A linear genome with distinct ends gives a path. A circular genome, or a linear one whose first and last (k−1)-mers coincide, gives a cycle.

## Suggested commit order
1. `euler.py`: Euler conditions + `check_euler`
2. Iterative Hierholzer + `path_to_sequence`
3. `unitigs` + ambiguity enumeration
4. `assembler.py` pipeline + `k_sweep`
5. `tests/test_part3_euler_assembly.py`
6. `demo.py`, `simulate.py`, `benchmark.py`
