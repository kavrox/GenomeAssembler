# Part 4 — Generalized Suffix Tree & Mutation Analysis

**Files:** `genome/suffix_tree.py`, `genome/mutations.py`
**Tests:** `tests/test_part4_gst_mutations.py`
**Guide sections:** Parts 12–16

## What this part does
It indexes the reference and the query together in one **generalized suffix tree** (GST) built
over `ref + '#' + query + '$'`. It finds their longest common substring and uses shared blocks as
**anchors**. Whatever lies between anchors is classified as a **substitution, insertion, deletion
or complex change**, with a 1-based position, HGVS-style notation and context.

## 1. Data structure justification
**Why a suffix tree?** Mutation analysis asks a *substring* question: which stretches do the two
sequences share, and where? Every substring is a prefix of some suffix. A suffix tree stores all
suffixes in compressed form, so:

| Question | Suffix tree | Alternative |
|---|---|---|
| Is P a substring? | O(m), independent of text length | naive scan O(n·m); KMP O(n+m) per query |
| All occurrences | O(m + occ) | rescan the text each time |
| Longest common substring of A and B | O(N) | dynamic programming O(|A|·|B|) |
| Longest repeat | O(N) | O(N²) brute force |

**Why *generalized*?** We compare two strings, and one tree answers questions about both. Each
string gets its **own unique sentinel** (`#`, `$`). Without them, a suffix could be a prefix of
another suffix and end in the middle of the tree with no leaf, and suffixes of A and B could merge.
With sentinels, every suffix ends at its own leaf labelled (string id, offset).

**Why compressed (not a trie)?** A suffix *trie* has O(n²) nodes. Path compression stores each
edge as a `(start, end)` pair of indices into the text, giving at most n leaves and fewer than n
internal nodes, so O(n) nodes.

**Honest comparison (guide Part 16):** for one pair of strings with one SNV, trimming the common
prefix and suffix is O(n) and needs almost no memory (our `classify_difference` does exactly that).
Suffix arrays and FM-indexes are far more memory-efficient, and real tools use them. The GST is
chosen because one index answers many questions (search, all occurrences, LCS, repeats), it scales
to several sequences, and it is the String-ADT structure the project is meant to demonstrate.

## 2. Implementation strategy
- Node = `children` dict keyed by first character (the 4 bases plus sentinels), edge label `(start, end)`, `leaf_info = (seq id, offset)`, `depth` = path-label length. `__slots__` keeps nodes small.
- **Two builders, same tree** (`test_naive_equals_ukkonen` checks this on 44 inputs):
  - **naive**: insert each suffix from the root, walk down, and split the edge at the first mismatch. Simple and easy to explain.
  - **Ukkonen**: online, one character at a time, with the *active point* (node, edge, length), *suffix links* and a *global end* for leaves (rule 1 grows every leaf for free). Rule 2 creates a leaf or splits; rule 3 means "already present, stop this phase".
- After building, one pass closes leaf ends, computes depths and labels leaves: `start = n − depth`.
- **LCS** is a post-order pass that computes, for every node, the leaf count and first offset per string. The answer is the deepest internal node with leaves from *both* strings. Ties go to the earliest position, so the result is deterministic. Internal nodes can never contain a sentinel, because a unique character can't occur twice.
- **Anchor-gap algorithm** (`anchor_alignment`): take the LCS as an anchor, then split into left and right sub-problems. It uses an **explicit stack, not recursion**, so very different sequences can't hit the recursion limit. When the best remaining shared block is shorter than `min_anchor` (3), the leftover pair (ref gap, query gap) is classified:
  - equal length → one SUBSTITUTION per differing base
  - empty ref gap → INSERTION (reported *before* that reference position)
  - empty query gap → DELETION
  - otherwise trim the common prefix and suffix, then INS/DEL, or COMPLEX (delins)
- `identity` = anchored fraction. Below 50% the report warns **LOW SIMILARITY**, because the variant list would not be meaningful.

## 3. Time complexity (with justification)
- **Naive build O(N²) worst case.** Inserting suffix i walks down up to (N − i) characters, and Σ(N − i) = N(N+1)/2. On random DNA it is about O(N log N), since walks only go as deep as the longest repeat (~log₄ N). Measured on `A^n`: ≈4× per doubling.
- **Ukkonen O(N)** for a constant alphabet: N phases; the total number of explicit extensions is ≤ N because `remaining` increases by 1 per phase and each extension decreases it; suffix links plus skip/count make each extension amortised O(1). Measured: ≈2× per doubling.
- **Space O(N) nodes**: at most N leaves and N − 1 internal nodes (each internal node has ≥ 2 children). N = |ref| + |query| + 2.
- **Search O(m + occ)**: walk m characters from the root, then collect the leaves below. The measured time per query barely changes as the text grows 16×.
- **LCS O(N)**: one post-order traversal.
- **Variant detection** builds one GST per segment. Each recursion level covers ≤ N characters in total, so it is O(N·d) with d = depth (≈ log of the number of anchors when splits are balanced), and O(N²) in the worst case. Gap classification is O(gap length).

## 4. Live demo (my part)
```bash
python main.py compare --ref ATGCGATAC --query ATGCAATAC          # SUBSTITUTION pos 5 G>A
python main.py compare --ref ATGCGATAC --query ATGCGGATAC         # INSERTION before pos 6 (g.5_6insG)
python main.py compare --ref ATGCGATAC --query ATGCATAC           # DELETION pos 5
python main.py compare --ref ATGCGATACCGTA --query ATGCAATACCGGA  # two SNVs
python main.py compare --ref AAAA --query TTTT                    # no anchor -> LOW SIMILARITY
python main.py tree BANANA                                        # print a suffix tree
python main.py search --text ATGCGATAC ATGCAATAC --pattern ATAC   # seq1:6, seq2:6 (1-based)
python main.py compare --ref data/sim/reference.fasta --query data/sim/query.fasta --no-alignment
python main.py benchmark                                          # naive O(n^2) vs Ukkonen O(n)
```

## 5. Questions I should be ready for
- *Suffix vs substring?* A suffix runs to the end; every substring is a prefix of some suffix.
- *Why sentinels?* So no suffix is a prefix of another, every suffix gets its own leaf, and leaves know which string they came from.
- *GST vs suffix tree?* A GST indexes several strings with unique sentinels and string-labelled leaves.
- *How do you find the LCS?* Take the deepest internal node whose subtree has leaves from both strings.
- *How are mutations found?* LCS anchors, recurse left and right, then classify the gaps.
- *Can it detect every mutation?* No. It handles clean SNVs and small indels; it can't find structural variants or rearrangements; variants closer together than `min_anchor` merge into one COMPLEX call; indels inside repeats have ambiguous positions.
- *Is a variant a disease mutation?* No. The tree knows nothing about biology. Say **SNV, not SNP** (SNP needs population frequency).
- *Why not just compare the strings?* For one SNV a linear scan is enough. The GST pays off for repeated queries, all-occurrence searches and LCS/repeat questions.
- *Naive vs Ukkonen?* Same tree. Naive is O(N²) worst case, Ukkonen O(N), and `benchmark` section 3 shows the difference.

## Suggested commit order
1. `Node` + naive GST construction + `contains` / `occurrences`
2. Post-order statistics: LCS, longest repeat, dump
3. Ukkonen construction + naive/Ukkonen equivalence test
4. `mutations.py`: gap classification, anchor alignment, variants, alignment rows
5. `MutationAnalyzer` / `MutationReport`
6. `tests/test_part4_gst_mutations.py`
