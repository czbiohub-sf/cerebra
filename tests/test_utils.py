import unittest

import pandas as pd
import vcfpy

from cerebra.utils import GenomePosition, GenomeIntervalTree

class GenomePositionTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.GPOS_A = GenomePosition("1", 0, 10)
        self.GPOS_B = GenomePosition("2", 5, 6)
        self.GPOS_C = GenomePosition("test", 100, 101)
        self.GPOS_D = GenomePosition("1", 5, 6)
        self.GPOS_E = GenomePosition("2", 0, 10)

        self.GPOS_A_STR = "1:1-10"
        self.GPOS_B_STR = "2:6-6"
        self.GPOS_C_STR = "test:101-101"

        self.GPOS_A_LEN = 10
        self.GPOS_B_LEN = 1

    def test_from_str(self):
        tests = [
            (self.GPOS_A_STR, self.GPOS_A),
            (self.GPOS_B_STR, self.GPOS_B),
            (self.GPOS_C_STR, self.GPOS_C),
        ]

        for test, expected in tests:
            self.assertEqual(expected, GenomePosition.from_str(test))

    def test_from_vcf_record(self):
        tests = [
            self.GPOS_A,
            self.GPOS_B,
            self.GPOS_C,
            self.GPOS_D,
            self.GPOS_E,
        ]

        tests = [
            (vcfpy.Record(
                CHROM=pos.chrom, POS=pos.start + 1, ID='.', REF='.' * len(pos),
                ALT=[], QUAL=0, FILTER='.', INFO={}
            ), pos) for pos in tests
        ]

        for test, expected in tests:
            self.assertEqual(expected, GenomePosition.from_vcf_record(test))

    def test_from_gtf_record(self):
        tests = [
            (pd.Series(["1", "unknown", "exon", "1", "10"]), self.GPOS_A),
            (pd.Series(["2", "unknown", "stop codon", "6", "6"]), self.GPOS_B),
            (pd.Series(["test", "unknown", "exon", "101", "101"]), self.GPOS_C),
        ]

        for test, expected in tests:
            self.assertEqual(expected, GenomePosition.from_gtf_record(test))

    def test_contains(self):
        positive_tests = [
            (self.GPOS_A, self.GPOS_D),
            (self.GPOS_E, self.GPOS_B),
            (self.GPOS_A, self.GPOS_A),
            (self.GPOS_B, self.GPOS_B),
        ]

        negative_tests = [
            (self.GPOS_A, self.GPOS_B),
            (self.GPOS_E, self.GPOS_D),
            (self.GPOS_D, self.GPOS_A),
            (self.GPOS_B, self.GPOS_E),
        ]

        for outer, inner in positive_tests:
            self.assertTrue(outer.contains(inner))

        for outer, inner in negative_tests:
            self.assertFalse(outer.contains(inner))

    def test_eq(self):
        positive_tests = [
            (self.GPOS_A, self.GPOS_A),
            (self.GPOS_B, self.GPOS_B),
            (self.GPOS_C, self.GPOS_C),
        ]

        negative_tests = [
            (self.GPOS_A, self.GPOS_B),
            (self.GPOS_B, self.GPOS_A),
            (self.GPOS_D, self.GPOS_A),
            (self.GPOS_E, self.GPOS_B),
            (self.GPOS_C, self.GPOS_B),
        ]

        for left, right in positive_tests:
            self.assertEqual(left, right)

        for left, right in negative_tests:
            self.assertNotEqual(left, right)

    def test_str(self):
        tests = [
            (self.GPOS_A, self.GPOS_A_STR),
            (self.GPOS_B, self.GPOS_B_STR),
        ]

        for test, expected in tests:
            self.assertEqual(expected, str(test))

    def test_len(self):
        tests = [
            (self.GPOS_A, self.GPOS_A_LEN),
            (self.GPOS_B, self.GPOS_B_LEN),
        ]

        for test, expected in tests:
            self.assertEqual(expected, len(test))

class GenomeIntervalTreeTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.tree_positions = [
            GenomePosition("1", 0, 10),
            GenomePosition("1", 20, 30),
            GenomePosition("1", 40, 41),
            GenomePosition("1", 50, 50),
            GenomePosition("2", 0, 10),
            GenomePosition("X", 100, 200),
        ]

        self.tree = GenomeIntervalTree(lambda row: row, self.tree_positions)

        # Multiple matches are specified because the first match can be any
        # matching interval. Similarly, multiple "best" matches are specified
        # because if two matches have the same Jaccard index, either can be
        # selected by the implementation.

        # (query, bests, matches)
        self.overlap_positive_tests = [
            (GenomePosition("1", 0, 10), self.tree_positions[0:1], self.tree_positions[0:1]),
            (GenomePosition("1", 9, 10), self.tree_positions[0:1], self.tree_positions[0:1]),
            (GenomePosition("2", -10, 20), self.tree_positions[4:5], self.tree_positions[4:5]),
            (GenomePosition("1", 40, 41), self.tree_positions[2:3], self.tree_positions[2:3]),
            (GenomePosition("1", 49, 51), self.tree_positions[3:4], self.tree_positions[3:4]),
            (GenomePosition("1", 0, 51), self.tree_positions[0:2], self.tree_positions[0:4]),
        ]

        self.overlap_negative_tests = [
            GenomePosition("1", 10, 19),
            GenomePosition("1", 19, 20),
            GenomePosition("2", 20, 30),
            GenomePosition("1", 20, 20),
            GenomePosition("1", 50, 50),
            GenomePosition("1", 49, 50),
            GenomePosition("1", 50, 51),
        ]

        # (query, bests, matches)
        self.containment_positive_tests = [
            (GenomePosition("1", 0, 10), self.tree_positions[0:1], self.tree_positions[0:1]),
            (GenomePosition("1", 0, 1), self.tree_positions[0:1], self.tree_positions[0:1]),
            (GenomePosition("1", 9, 10), self.tree_positions[0:1], self.tree_positions[0:1]),
            (GenomePosition("1", 2, 8), self.tree_positions[0:1], self.tree_positions[0:1]),
            (GenomePosition("1", 40, 41), self.tree_positions[2:3], self.tree_positions[2:3]),
        ]

        self.containment_negative_tests = [
            GenomePosition("1", 0, 30),
            GenomePosition("1", 0, 11),
            GenomePosition("2", 22, 28),
            GenomePosition("1", 50, 50),
            GenomePosition("1", 40, 40),
            GenomePosition("1", 10, 11),
        ]

    def test_compute_jaccard_index(self):
        tests = [
            ((GenomePosition("1", 0, 10), GenomePosition("1", 5, 15)), 5 / 15),
            ((GenomePosition("1", 100, 0), GenomePosition("1", 100, 0)), 0),
            ((GenomePosition("1", 0, 1), GenomePosition("2", 0, 1)), 0),
            ((GenomePosition("1", 0, 1), GenomePosition("1", 0, 1)), 1),
        ]

        for test, expected in tests:
            self.assertAlmostEqual(expected, self.tree._compute_jaccard_index(*test))

    def test_has_overlap(self):
        for test, _, _ in self.overlap_positive_tests:
            self.assertTrue(self.tree.has_overlap(test))

        for test in self.overlap_negative_tests:
            self.assertFalse(self.tree.has_overlap(test))

    def test_overlap_finding(self):
        for test, expected_bests, expected_matches in self.overlap_positive_tests:
            first = self.tree.get_first_overlap(test)
            self.assertIsNotNone(first)
            self.assertIn(first, expected_matches)

            best = self.tree.get_best_overlap(test)
            self.assertIn(best, expected_bests)

            matches = self.tree.get_all_overlaps(test)
            self.assertEqual(expected_matches, matches)

        for test in self.overlap_negative_tests:
            first = self.tree.get_first_overlap(test)
            self.assertIsNone(first)

            matches = self.tree.get_all_overlaps(test)
            self.assertEqual([], matches)

    def test_containment_finding(self):
        for test, expected_bests, expected_matches in self.containment_positive_tests:
            first = self.tree.get_first_containment(test)
            self.assertIsNotNone(first)
            self.assertIn(first, expected_matches)

            best = self.tree.get_best_overlap(test)
            self.assertIn(best, expected_bests)

            matches = self.tree.get_all_containments(test)
            self.assertEqual(expected_matches, matches)

        for test in self.containment_negative_tests:
            first = self.tree.get_first_containment(test)
            self.assertIsNone(first)

            matches = self.tree.get_all_containments(test)
            self.assertEqual([], matches)


if __name__ == "__main__":
    unittest.main()
