"""Part 3 tests - Euler conditions, Hierholzer, unitigs, assembly pipeline."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genome.assembler import assemble, k_sweep
from genome.euler import (all_euler_sequences, check_euler, euler_endpoints, hierholzer,
                          path_to_sequence, unitigs)
from genome.graph import DeBruijnGraph
from genome.simulate import random_genome, simulate_reads

PART18_READS = ["ATGC", "GCGAT", "GATAC"]
ERROR_READS = ["ATGCGATAC"] * 3 + ["ATGCGTTAC"]


class TestEuler(unittest.TestCase):
    def test_endpoints_part18(self):
        g = assemble(PART18_READS, 3).graph
        self.assertEqual(euler_endpoints(g), ("AT", "AC", True))
        check = check_euler(g)
        self.assertTrue(check.ok)
        self.assertEqual(check.kind, "path")

    def test_hierholzer_trace(self):                               # guide Part 8.2
        g = assemble(PART18_READS, 3).graph
        path = hierholzer(g)
        self.assertEqual(path, ["AT", "TG", "GC", "CG", "GA", "AT", "TA", "AC"])
        self.assertEqual(len(path), g.num_edges() + 1)
        self.assertEqual(path_to_sequence(path), "ATGCGATAC")

    def test_cycle(self):                                          # Part 35 test 1
        r = assemble(["ATGCGAT"], 3)
        self.assertEqual(r.method, "eulerian_cycle")
        self.assertEqual(r.sequence, "ATGCGAT")

    def test_degree_violation(self):
        g = assemble(["ATGC", "CCGGA"], 3).graph
        self.assertEqual(euler_endpoints(g), (None, None, False))
        self.assertIsNone(hierholzer(g))
        self.assertIn("degree condition", check_euler(g).reason)

    def test_balanced_but_disconnected(self):
        g = DeBruijnGraph(3)
        for a, b in (("AC", "CA"), ("CA", "AC"), ("GT", "TG"), ("TG", "GT")):
            g.add_edge(a, b)
        self.assertTrue(euler_endpoints(g)[2])        # degrees balance...
        self.assertIsNone(hierholzer(g))              # ...but one walk cannot cover both
        self.assertIn("disconnected", check_euler(g).reason)

    def test_empty_graph(self):
        self.assertEqual(check_euler(DeBruijnGraph(3)).kind, "empty")
        self.assertEqual(path_to_sequence([]), "")

    def test_unitigs_include_cycles(self):
        g = DeBruijnGraph.from_sequence("ATGCGAT", 3, multigraph=False)
        self.assertEqual(unitigs(g), ["ATGCGAT"])


class TestAssembly(unittest.TestCase):
    def test_part18(self):                                         # Part 35 test 2
        r = assemble(PART18_READS, 3)
        self.assertEqual(r.method, "eulerian_path")
        self.assertEqual(r.sequence, "ATGCGATAC")
        self.assertIs(r.ambiguous, False)

    def test_duplicated_reads(self):                               # Part 35 test 3
        self.assertEqual(assemble(PART18_READS * 2, 3).sequence, "ATGCGATAC")

    def test_two_components(self):                                 # Part 35 test 8
        r = assemble(["ATGC", "CCGGA"], 3)
        self.assertEqual(r.method, "unitigs")
        self.assertEqual(r.contigs, ["ATGC", "CCGGA"])

    def test_error_bubble(self):                                   # Part 35 test 9
        r = assemble(ERROR_READS, 3)
        self.assertIsNone(r.sequence)
        self.assertEqual(r.contigs, ["ATA", "ATGCG", "CGAT", "CGTTA", "TAC"])
        self.assertEqual(assemble(ERROR_READS, 3, min_coverage=2).sequence, "ATGCGATAC")
        cleaned = assemble(ERROR_READS, 3, clean=True)
        self.assertEqual(cleaned.sequence, "ATGCGATAC")
        self.assertEqual(cleaned.cleaning["bubble_edges_removed"], 3)

    def test_k_too_large(self):                                    # Part 35 test 10
        r = assemble(PART18_READS, 4)
        self.assertEqual((r.graph.num_vertices(), r.graph.num_edges()), (7, 5))
        self.assertEqual(r.contigs, ["ATGC", "GCGATAC"])

    def test_k_too_small(self):                                    # guide Part 10.2
        r = assemble(PART18_READS, 2)
        self.assertEqual(r.sequence, "TACGATGC")
        self.assertTrue(r.ambiguous)

    def test_repeat_ambiguity(self):                               # Part 35 tests 11 + 12
        g3 = DeBruijnGraph.from_sequence("TTGTGAGTA", 3, multigraph=True)
        self.assertEqual(all_euler_sequences(g3), {"TTGTGAGTA", "TTGAGTGTA"})
        g4 = DeBruijnGraph.from_sequence("TTGTGAGTA", 4, multigraph=True)
        self.assertEqual(all_euler_sequences(g4), {"TTGTGAGTA"})

    def test_k_table_part10(self):                                 # guide Part 10.1
        expected = {2: (3, 8, 15), 3: (6, 7, 2), 4: (7, 6, 1), 5: (6, 5, 1), 9: (2, 1, 1)}
        for k, (v, e, n) in expected.items():
            g = DeBruijnGraph.from_sequence("TTGTGAGTA", k, multigraph=True)
            self.assertEqual((g.num_vertices(), g.num_edges(), len(all_euler_sequences(g))), (v, e, n))

    def test_compeau_example(self):                                # guide Part 37
        r = assemble(["ATGGCGTGCA"], 3)
        self.assertEqual(r.sequence, "ATGCGTGGCA")
        self.assertTrue(r.ambiguous)
        self.assertIn("ATGGCGTGCA", r.alternatives)

    def test_distinct_edges_collapse_repeats(self):                # guide Part 34
        self.assertEqual(assemble(["ATGCGAT"], 2).sequence, "ATGCGA")
        self.assertEqual(assemble(["ATGCGAT"], 2, multigraph=True).sequence, "ATGCGAT")

    def test_k_longer_than_reads(self):                            # Part 35 test 13
        r = assemble(["ATG"], 5)
        self.assertEqual(r.method, "empty")
        self.assertEqual(r.contigs, [])
        self.assertIsNone(r.sequence)

    def test_invalid_input(self):                                  # Part 35 test 14
        with self.assertRaises(ValueError):
            assemble(["ATXG"], 3)
        with self.assertRaises(ValueError):
            assemble(PART18_READS, 1)
        with self.assertRaises(ValueError):
            assemble([], 3)

    def test_k_sweep(self):
        rows = {r["k"]: r for r in k_sweep(PART18_READS, [2, 3, 4])}
        self.assertEqual(rows[3]["sequence"], "ATGCGATAC")
        self.assertEqual(rows[4]["method"], "unitigs")
        self.assertEqual(rows[2]["sequence"], "TACGATGC")


class TestSimulatedGenomes(unittest.TestCase):
    def test_error_free_reads_reconstruct_genome(self):
        for seed in range(3):
            genome = random_genome(1500, seed)
            reads = simulate_reads(genome, 80, 8, seed=seed)
            r = assemble(reads, 25, check_ambiguity=False)
            self.assertEqual(r.sequence, genome)

    def test_noisy_reads_with_filter(self):
        genome = random_genome(1500, 11)
        reads = simulate_reads(genome, 100, 30, error_rate=0.005, seed=11)
        r = assemble(reads, 25, min_coverage=3, check_ambiguity=False)
        self.assertTrue(r.success)
        self.assertIn(r.sequence, genome)               # exact, only low-coverage ends trimmed
        self.assertGreater(len(r.sequence), 0.97 * len(genome))


if __name__ == "__main__":
    unittest.main()
