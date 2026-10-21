"""Part 2 tests - Graph ADT / De Bruijn graph construction and cleaning."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genome.graph import DeBruijnGraph
from genome.kmers import count_kmers

PART18_READS = ["ATGC", "GCGAT", "GATAC"]
ERROR_READS = ["ATGCGATAC"] * 3 + ["ATGCGTTAC"]


def graph_of(reads, k, **kw):
    return DeBruijnGraph.from_kmer_counts(count_kmers(reads, k), k, **kw)


class TestConstruction(unittest.TestCase):
    def test_single_read_cycle(self):                              # Part 35 test 1
        g = graph_of(["ATGCGAT"], 3)
        self.assertEqual((g.num_vertices(), g.num_edges()), (5, 5))
        self.assertTrue(all(d == (1, 1) for d in g.degree_report().values()))

    def test_part18_graph(self):                                   # Part 35 test 2
        g = graph_of(PART18_READS, 3)
        self.assertEqual((g.num_vertices(), g.num_edges()), (7, 7))
        self.assertEqual(g.degree_report()["AT"], (1, 2))
        self.assertEqual(g.degree_report()["AC"], (1, 0))
        self.assertEqual(sorted(g.neighbors("AT")), ["TA", "TG"])
        self.assertEqual(g.branching_vertices(), ["AT"])
        gat = next(e for e in g.edges() if e.kmer == "GAT")
        self.assertEqual((gat.src, gat.dst, gat.coverage), ("GA", "AT", 2))

    def test_vertices_are_k_minus_1_mers(self):
        g = graph_of(["ATGCGATACCA"], 5)
        self.assertTrue(all(len(v) == 4 for v in g.vertices))
        self.assertTrue(all(e.kmer == e.src + e.dst[-1] for e in g.edges()))

    def test_out_degree_at_most_4(self):
        g = graph_of(["ACGTACGGACTTAGCAAGCTAGCATCGATCGACTGA"], 3)
        self.assertTrue(all(g.out_degree(v) <= 4 for v in g.vertices))

    def test_adt_operations(self):
        g = DeBruijnGraph(3)
        self.assertTrue(g.is_empty())
        e = g.add_edge("AT", "TG")
        self.assertEqual(e.kmer, "ATG")
        self.assertTrue(g.has_edge("AT", "TG"))
        self.assertFalse(g.has_edge("TG", "AT"))
        self.assertEqual((g.in_degree("TG"), g.out_degree("AT")), (1, 1))
        g.remove_edge(e)
        self.assertEqual(g.num_edges(), 0)
        self.assertEqual(g.remove_isolated_vertices(), 2)
        self.assertTrue(g.is_empty())

    def test_remove_vertex(self):
        g = graph_of(PART18_READS, 3)
        g.remove_vertex("AT")
        self.assertFalse(g.has_vertex("AT"))
        self.assertEqual(g.num_edges(), 4)       # ATA, ATG, GAT gone

    def test_edge_validation(self):
        g = DeBruijnGraph(3)
        with self.assertRaises(ValueError):
            g.add_edge("AT", "GC")               # no k-2 overlap
        with self.assertRaises(ValueError):
            g.add_edge("ATG", "TGC")             # wrong vertex length
        with self.assertRaises(ValueError):
            DeBruijnGraph(1)

    def test_multigraph_keeps_repeats(self):                       # Part 35 test 11
        g = DeBruijnGraph.from_sequence("TTGTGAGTA", 3, multigraph=True)
        self.assertEqual(g.degree_report()["TG"], (2, 2))
        self.assertEqual(g.degree_report()["GT"], (2, 2))
        collapsed = DeBruijnGraph.from_sequence("ATGCGAT", 2, multigraph=False)
        self.assertEqual(collapsed.num_edges(), 5)                  # AT counted once

    def test_components(self):                                     # Part 35 test 8
        g = graph_of(["ATGC", "CCGGA"], 3)
        self.assertEqual(len(g.weakly_connected_components()), 2)
        self.assertEqual(g.summary()["components"], 2)

    def test_copy_is_independent(self):
        g = graph_of(PART18_READS, 3)
        h = g.copy()
        h.remove_vertex("AT")
        self.assertEqual(g.num_edges(), 7)


class TestCleaning(unittest.TestCase):
    def test_coverage_filter(self):                                # Part 35 test 9
        g = graph_of(ERROR_READS, 3)
        self.assertEqual(g.num_edges(), 10)
        removed = g.filter_low_coverage(2)
        self.assertEqual(removed, 3)
        self.assertEqual(sorted(e.kmer for e in g.edges()),
                         ["ATA", "ATG", "CGA", "GAT", "GCG", "TAC", "TGC"])
        filtered = graph_of(ERROR_READS, 3, min_coverage=2)
        self.assertEqual(filtered.num_edges(), 7)

    def test_bubble_popping(self):
        g = graph_of(ERROR_READS, 3)
        self.assertEqual(g.pop_bubbles(), 3)
        self.assertFalse(any(e.kmer in ("CGT", "GTT", "TTA") for e in g.edges()))
        self.assertFalse(g.has_vertex("GT"))

    def test_tip_removal(self):
        g = graph_of(["ATGCGATAC"] * 3 + ["ATGCGATAG"], 3)       # error in the last base
        self.assertTrue(g.has_vertex("AG"))
        self.assertEqual(g.remove_tips(), 1)
        self.assertFalse(g.has_vertex("AG"))

    def test_cleaning_keeps_real_low_coverage_structure(self):
        g = graph_of(PART18_READS, 3)            # every k-mer seen once: nothing is an error
        self.assertEqual(g.remove_tips(), 0)
        self.assertEqual(g.pop_bubbles(), 0)
        self.assertEqual(g.num_edges(), 7)

    def test_balanced_bubble_is_kept(self):
        # two alleles with equal coverage = a real variant, not an error
        g = graph_of(["ATGCGATAC"] * 2 + ["ATGCGTTAC"] * 2, 3)
        self.assertEqual(g.pop_bubbles(), 0)


class TestExport(unittest.TestCase):
    def test_dot(self):
        g = graph_of(PART18_READS, 3)
        path = ["AT", "TG", "GC", "CG", "GA", "AT", "TA", "AC"]
        dot = g.to_dot(path)
        self.assertIn('"AT" -> "TG" [label="ATG x1 [#1]"]', dot)
        self.assertIn('"TA" -> "AC" [label="TAC x1 [#7]"]', dot)
        self.assertTrue(dot.startswith("digraph"))


if __name__ == "__main__":
    unittest.main()
