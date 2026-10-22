"""Part 4 tests - Generalized suffix tree and mutation analysis."""

import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genome.mutations import (MutationAnalyzer, Variant, anchor_alignment, classify_difference,
                              classify_gap, find_variants)
from genome.simulate import mutate, random_genome
from genome.suffix_tree import GeneralizedSuffixTree, SuffixTree

ALGOS = ("naive", "ukkonen")


def calls(ref, qry, **kw):
    return [v.as_tuple() for v in find_variants(ref, qry, **kw)]


class TestSuffixTree(unittest.TestCase):
    def test_banana(self):                                         # Part 35 test 16
        for algo in ALGOS:
            t = SuffixTree("BANANA", algo)
            self.assertEqual(t.occurrences("ANA"), [(0, 1), (0, 3)])
            self.assertEqual(t.count_leaves(), 7)
            self.assertEqual(t.longest_repeated_substring(), "ANA")
            self.assertTrue(t.contains("NAN"))
            self.assertFalse(t.contains("NAB"))

    def test_part18_gst(self):                                     # guide Part 18.6 / demo
        for algo in ALGOS:
            t = GeneralizedSuffixTree(["ATGCGATAC", "ATGCAATAC"], algo)
            self.assertEqual(len(t.text), 20)
            self.assertEqual(t.count_nodes(), 32)
            self.assertEqual(t.lcs_with_positions(), ("ATGC", {0: 0, 1: 0}))
            self.assertEqual(t.occurrences("ATAC"), [(0, 5), (1, 5)])

    def test_lcs_part13(self):                                     # Part 35 test 17
        for algo in ALGOS:
            t = GeneralizedSuffixTree(["ATCGATCG", "ATGGATCG"], algo)
            self.assertEqual(t.longest_common_substring(), "GATCG")
            self.assertEqual(t.occurrences("ATCG"), [(0, 0), (0, 4), (1, 4)])

    def test_naive_equals_ukkonen(self):
        rng = random.Random(3)
        cases = [["A" * 30, "A" * 25], ["ABABABAB"], ["MISSISSIPPI", "PIPPI"],
                 ["ACGT" * 10, "ACGA" * 8, "TTTT"]]
        for _ in range(40):
            cases.append(["".join(rng.choice("ACG") for _ in range(rng.randint(1, 40)))
                          for _ in range(rng.randint(1, 3))])
        for seqs in cases:
            naive = GeneralizedSuffixTree(seqs, "naive")
            ukk = GeneralizedSuffixTree(seqs, "ukkonen")
            self.assertEqual(naive.signature(), ukk.signature(), seqs)
            total = sum(len(s) + 1 for s in seqs)
            self.assertEqual(naive.count_leaves(), total)

    def test_occurrences_match_brute_force(self):
        rng = random.Random(5)
        text = "".join(rng.choice("ACGT") for _ in range(300))
        tree = GeneralizedSuffixTree([text], "ukkonen")
        for _ in range(100):
            p = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 4)))
            brute = [(0, i) for i in range(len(text)) if text.startswith(p, i)]
            self.assertEqual(tree.occurrences(p), brute)

    def test_lcs_matches_brute_force(self):
        rng = random.Random(9)
        for _ in range(30):
            a = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 25)))
            b = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 25)))
            best = max((len(a[i:j]) for i in range(len(a)) for j in range(i + 1, len(a) + 1)
                        if a[i:j] in b), default=0)
            self.assertEqual(len(GeneralizedSuffixTree([a, b]).longest_common_substring()), best)

    def test_input_checks(self):
        with self.assertRaises(ValueError):
            GeneralizedSuffixTree(["AC#GT"])
        with self.assertRaises(ValueError):
            GeneralizedSuffixTree(["ACGT"], "quadratic")
        with self.assertRaises(ValueError):
            GeneralizedSuffixTree([])
        t = GeneralizedSuffixTree(["ACGT", "TTT"])
        self.assertEqual(t.occurrences("#"), [])
        self.assertEqual(t.occurrences("GG"), [])

    def test_dump(self):
        dump = SuffixTree("BANANA").dump()
        self.assertIn("+-- BANANA#   <- leaf (seq1, pos 0)", dump)
        self.assertTrue(dump.startswith("root"))


class TestMutations(unittest.TestCase):
    def test_substitution(self):                                   # Part 35 test 2
        self.assertEqual(calls("ATGCGATAC", "ATGCAATAC"), [("SUBSTITUTION", 5, "G", "A")])
        self.assertEqual(calls("ATGCTAGCTA", "ATGCCAGCTA"), [("SUBSTITUTION", 5, "T", "C")])

    def test_insertion(self):                                      # Part 35 test 4
        self.assertEqual(calls("ATGCGATAC", "ATGCGGATAC"), [("INSERTION", 6, "", "G")])

    def test_deletion(self):                                       # Part 35 test 5
        self.assertEqual(calls("ATGCGATAC", "ATGCATAC"), [("DELETION", 5, "G", "")])

    def test_two_snvs(self):                                       # Part 35 test 6
        self.assertEqual(calls("ATGCGATACCGTA", "ATGCAATACCGGA"),
                         [("SUBSTITUTION", 5, "G", "A"), ("SUBSTITUTION", 12, "T", "G")])

    def test_identical(self):                                      # Part 35 test 7
        self.assertEqual(calls("ATGCGATAC", "ATGCGATAC"), [])

    def test_no_anchor(self):                                      # Part 35 test 15
        self.assertEqual(len(calls("AAAA", "TTTT")), 4)
        report = MutationAnalyzer("AAAA", "TTTT").report()
        self.assertTrue(report.low_similarity)
        self.assertEqual(report.identity, 0)

    def test_part13_substitution(self):                            # Part 35 test 17
        self.assertEqual(calls("ATCGATCG", "ATGGATCG"), [("SUBSTITUTION", 3, "C", "G")])

    def test_part15_example(self):
        self.assertEqual(calls("ACGTACGT", "ACGTTCGT"), [("SUBSTITUTION", 5, "A", "T")])

    def test_ukkonen_gives_same_calls(self):
        for ref, qry in (("ATGCGATAC", "ATGCGGATAC"), ("ATGCGATACCGTA", "ATGCAATACCGGA")):
            self.assertEqual(calls(ref, qry, algorithm="ukkonen"), calls(ref, qry))

    def test_classify_gap(self):
        self.assertEqual(classify_gap("", "GG", 4), [Variant("INSERTION", 4, "", "GG")])
        self.assertEqual(classify_gap("TT", "", 4), [Variant("DELETION", 4, "TT", "")])
        self.assertEqual([v.as_tuple() for v in classify_gap("ACG", "AG", 1)], [("DELETION", 2, "C", "")])
        self.assertEqual([v.as_tuple() for v in classify_gap("AT", "GCC", 7)], [("COMPLEX", 7, "AT", "GCC")])
        self.assertEqual(classify_gap("A", "A", 1), [])

    def test_single_edit_shortcut(self):
        self.assertEqual([v.as_tuple() for v in classify_difference("ATGCGATAC", "ATGCATAC")],
                         [("DELETION", 5, "G", "")])

    def test_hgvs_and_context(self):
        vs = find_variants("ATGCGATAC", "ATGCAATAC")
        self.assertEqual((vs[0].hgvs, vs[0].context), ("g.5G>A", "ATGC[G>A]ATAC"))
        self.assertEqual(find_variants("ATGCGATAC", "ATGCGGATAC")[0].hgvs, "g.5_6insG")
        self.assertEqual(find_variants("ATGCGATAC", "ATGCATAC")[0].hgvs, "g.5del")
        self.assertEqual(Variant("DELETION", 3, "GAT", "").hgvs, "g.3_5del")

    def test_alignment(self):
        rep = MutationAnalyzer("ATGCGATAC", "ATGCGGATAC").report()
        self.assertEqual(rep.alignment(), ("ATGCG-ATAC", "||||| ||||", "ATGCGGATAC"))
        rep = MutationAnalyzer("ATGCGATAC", "ATGCAATAC").report()
        self.assertEqual(rep.alignment()[1], "||||x||||")
        self.assertAlmostEqual(rep.identity, 8 / 9)

    def test_anchors_are_ordered_and_disjoint(self):
        anchors, gaps = anchor_alignment("ATGCGATACCGTA", "ATGCAATACCGGA")
        self.assertEqual([(a.ref_start, a.query_start, a.length) for a in anchors],
                         [(0, 0, 4), (5, 5, 6)])
        self.assertEqual(gaps, [(4, 5, 4, 5), (11, 13, 11, 13)])
        for a, b in zip(anchors, anchors[1:]):
            self.assertLessEqual(a.ref_end, b.ref_start)
            self.assertLessEqual(a.query_end, b.query_start)

    def test_planted_mutations_are_recovered(self):
        for seed in range(6):
            ref = random_genome(1200, seed)
            qry, truth = mutate(ref, n_snv=3, n_ins=2, n_del=2, seed=seed)
            self.assertEqual(calls(ref, qry), truth, f"seed {seed}")

    def test_validation(self):
        with self.assertRaises(ValueError):
            MutationAnalyzer("ATGN", "ATGC")
        with self.assertRaises(ValueError):
            MutationAnalyzer("ATGC", "ATGC", min_anchor=0)


if __name__ == "__main__":
    unittest.main()
