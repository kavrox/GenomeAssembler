# In-Memory Genomic Sequence Assembler & Mutation Analyzer

**Graph ADT + String ADT · De Bruijn Graph · Eulerian path (Hierholzer) · Generalized Suffix Tree (naive + Ukkonen)**

Team: Anush Sivakumar · Arjun Padmakumar · Shris · Theajaa

The program takes short overlapping DNA reads and rebuilds the original sequence
(de novo assembly). It then compares a sequence with a reference and reports
candidate substitutions, insertions and deletions (mutation analysis).

```
          GRAPH ADT PIPELINE                                 STRING ADT PIPELINE
reads -> k-mers -> De Bruijn graph -> Eulerian path -> CONTIG -> GST(ref # query $) -> anchors/gaps -> variants
         (hash map)  (adjacency list)    (Hierholzer)   (string)   (suffix tree)        (LCS)            (report)
```

The two halves share only a plain string (the assembled contig), so either one
can be swapped out without touching the other.

---

## 1. Running it

Python 3.8+ only, no third-party packages. (On the project laptop, `python` isn't on PATH, so use
`C:\msys64\ucrt64\bin\python3.exe` wherever this README says `python`.)

```bash
python main.py demo                 # the 8-step walk-through (k=3), < 1 s
python main.py demo 4               # k too large  -> graph fragments, no Euler path
python main.py demo 2               # k too small  -> valid Euler path, WRONG sequence
python main.py gui                  # graphical interface (Tkinter)
python -m unittest -v               # 96 tests
```

Any `--reads`, `--ref` or `--query` value can be a **FASTA / FASTQ / plain-text file** or a **literal DNA
string**, so you can type faculty test cases straight into the command line:

```bash
# full pipeline: assemble, check against the reference, call variants
python main.py pipeline --reads ATGC GCGAT GATAC -k 3 --ref ATGCGATAC --query ATGCAATAC

# assembly only (sequencing-error example: filter or clean the graph)
python main.py assemble --reads data/error_reads.fasta -k 3                  # bubble -> unitigs
python main.py assemble --reads data/error_reads.fasta -k 3 --min-coverage 2 # fixed by filter
python main.py assemble --reads data/error_reads.fasta -k 3 --clean          # fixed by bubble popping
python main.py assemble --reads data/reads.fastq -k 3 --min-quality 20       # FASTQ quality filter

# mutation analysis only
python main.py compare --ref ATGCGATAC --query ATGCGGATAC                    # insertion
python main.py compare --ref data/reference.fasta --query data/query_deletion.fasta --algorithm ukkonen

# experiments for the viva
python main.py ksweep --reads data/reads.fasta --k 2-5                        # choosing k
python main.py search --text ATGCGATAC ATGCAATAC --pattern ATAC GAT           # GST queries
python main.py tree BANANA                                                     # print a suffix tree
python main.py benchmark                                                       # measured Big-O

# a bigger, realistic run: 2 kb genome, 300 reads, 8 planted mutations
python main.py pipeline --reads data/sim/reads.fasta -k 21 --ref data/sim/reference.fasta \
                        --query data/sim/query.fasta --quiet --profile --out out
python main.py simulate --length 5000 --error-rate 0.005 --coverage 30 --out data/sim5k   # new data
```

`--out DIR` writes `contigs.fasta`, `variants.tsv`, `graph.dot` (render with Graphviz:
`dot -Tpng out/graph.dot -o graph.png`) and the full text report.

| Option | Meaning |
|---|---|
| `-k` | k-mer length (vertices are (k-1)-mers, edges are k-mers) |
| `--min-coverage N` | drop k-mers seen fewer than N times (sequencing-error filter) |
| `--clean` | remove tips and pop bubbles in the graph |
| `--multigraph` | one edge per k-mer *occurrence* (only for a single known sequence) |
| `--min-quality Q` | drop FASTQ reads whose mean Phred quality is below Q |
| `--min-anchor N` | shortest shared block trusted as an anchor (default 3) |
| `--algorithm naive\|ukkonen` | suffix-tree construction algorithm |
| `--profile` | stage timings + peak memory |

---

## 2. Project layout and the 4-part split

Each member owns one part: its modules, its tests and its viva sheet in `docs/`.
Everyone commits from **their own GitHub account**.

| Part | Owner | Files | Data structures | Viva sheet |
|---|---|---|---|---|
| **1. String ADT, input, k-mers, reporting** | _(fill in)_ | `genome/sequence.py`, `genome/io.py`, `genome/kmers.py`, `genome/report.py`, `main.py`, `data/` | String ADT, hash map (Counter) | [docs/PART1_strings_kmers.md](docs/PART1_strings_kmers.md) |
| **2. De Bruijn graph** | _(fill in)_ | `genome/graph.py`, `gui.py` | Graph ADT, adjacency list (hash map of dynamic arrays) | [docs/PART2_debruijn_graph.md](docs/PART2_debruijn_graph.md) |
| **3. Eulerian assembly** | _(fill in)_ | `genome/euler.py`, `genome/assembler.py`, `genome/simulate.py`, `genome/benchmark.py`, `demo.py` | stack, edge-pointer array | [docs/PART3_euler_assembly.md](docs/PART3_euler_assembly.md) |
| **4. GST + mutation analysis** | _(fill in)_ | `genome/suffix_tree.py`, `genome/mutations.py` | generalized suffix tree | [docs/PART4_gst_mutations.md](docs/PART4_gst_mutations.md) |

Tests: `tests/test_part1_strings_kmers.py`, `tests/test_part2_graph.py`,
`tests/test_part3_euler_assembly.py`, `tests/test_part4_gst_mutations.py`, plus the shared
`tests/test_cli.py` (end-to-end).

### Interfaces between parts (do not change without telling the team)

```python
# Part 1
validate_dna(seq) -> str                                   # raises ValueError
count_kmers(reads, k) -> Counter                           # {kmer: coverage}
# Part 2
DeBruijnGraph.from_kmer_counts(freq, k, min_coverage=1) -> DeBruijnGraph
graph.out[v] / graph.inc[v] -> list[Edge(src, dst, kmer, coverage)]
# Part 3
assemble(reads, k, min_coverage=1, clean=False) -> AssemblyResult   # .sequence, .contigs, .path
# Part 4
MutationAnalyzer(ref, query).report() -> MutationReport             # .variants, .lcs, .alignment()
find_variants(ref, query) -> list[Variant(type, position, ref, alt)]
```

Dependencies: Part 1 ← Part 2 ← Part 3, while Part 4 depends only on `validate_dna`.
Commit and merge in the order 1 → 2 → 3 → 4, then `main.py` / `gui.py` / `test_cli.py`.

### GitHub workflow

1. The repo owner creates the repository, adds the other three as collaborators and pushes
   `README.md`, `.gitignore`, `genome/__init__.py` and `tests/__init__.py`.
2. Each member makes a branch (`part1-strings-kmers`, `part2-graph`, `part3-euler`,
   `part4-gst`), commits their own files in small, working steps and opens a pull request.
3. Before the first commit on any laptop, check it is *your* identity:
   `git config user.name` and `git config user.email`.
4. After every merge, run `python -m unittest` on `main` before the next merge.

---

## 3. Design summary (for the evaluation rubric)

| Component | Data structure / algorithm | Time | Space |
|---|---|---|---|
| DNA sequence | String ADT (array-backed `str`) | index O(1), slice O(m) | O(n) |
| k-mer generation | sliding window | O(L·k) per read (L−k+1 windows) | O(L·k) |
| k-mer counting | hash map (`Counter`) | O(K·k) expected | O(D·k) |
| Assembly graph | De Bruijn graph, adjacency lists in hash maps | build O(D·k) | O(V+E), E ≤ 4V |
| Error removal | coverage filter / tip removal / bubble popping | O(E) / O(V·k) / bounded local BFS per branch vertex | O(V) |
| Contig assembly | Euler check + iterative Hierholzer (stack) | O(V+E) | O(V+E) |
| Fallback | unitig compaction | O(V+E) | O(V+E) |
| Sequence index | generalized suffix tree, naive / Ukkonen | O(N²) worst / O(N) | O(N) nodes |
| Pattern search | walk from the root | O(m + occ) | O(1) |
| LCS | deepest node with leaves from both strings | O(N) | O(N) |
| Mutation calling | LCS anchors (one GST per segment) + gap classification | O(N·d), d = anchor recursion depth (≈ log #anchors); O(N²) worst | O(N) |

K = total k-mers, D = distinct k-mers, N = |ref| + |query| + 2.

**Measured** (`python main.py benchmark`, this laptop): assembly stages and the Ukkonen build grow
≈2× when the input doubles (linear). On the worst-case input `A^n` the naive suffix-tree build
grows ≈4× (quadratic) while Ukkonen stays ≈2×. Search time per query stays almost flat as
the text grows 16× (it depends on the pattern length, not the text).

### Honest limitations
Reads from one strand only. Error handling is coverage filtering plus tip/bubble removal; there
is no full error correction. Repeat ambiguity is detected and reported but not resolved. The
anchor-gap variant caller is a heuristic, not an aligner: no structural variants, and indels
inside repeats have ambiguous positions. A "variant" here is a **candidate sequence
difference**, not a clinical finding. This is a teaching-scale tool (thousands of bases),
not a genome-scale assembler.

Background reading: `context pdf.pdf` (study guide, Parts 1–39) and the presentation deck.
