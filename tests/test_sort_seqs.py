#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Test suite for the sort_seqs.py module."""

from cogent.util.unit_test import TestCase, main
from nested_reference_otus.sort_seqs import compute_sequence_stats

class SortSeqsTests(TestCase):
    """Tests for the sort_seqs.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.fasta1 = [">1", "AGGT", ">2", "AGGA", ">3", "AGGC"]
        self.tax_map1 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tG1\tA;B;C\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;B;Z\tfoo\n"]
        self.tax_map2 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tG1\tA;B;C\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;B;Z;Z;F\tfoo\n"]
        # Taxonomy string with no semicolons (only one level).
        self.tax_map3 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tG1\tA;B;C\tfoo\n",
                "2\tG2\tA\tfoo\n",
                "3\tG3\tA;B;Z;Z;F\tfoo\n"]
        # Missing column in one of the rows.
        self.tax_map_invalid1 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tA;B;C\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;B;Z;Z;F\tfoo\n"]
        # Missing header.
        self.tax_map_invalid2 = [
                "1\tG1\tA;B;C\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;B;Z;Z;F\tfoo\n"]

    def test_compute_sequence_stats_unequal_depths(self):
        """Test computing seq stats on unequal taxonomic depths."""
        exp = {'1': [3, 4, 'AGGT'], '3': [2, 4, 'AGGC'], '2': [3, 4, 'AGGA']}
        obs = compute_sequence_stats(self.fasta1, self.tax_map1, ['Z'])
        self.assertEqual(obs, exp)

    def test_compute_sequence_stats_single_depth(self):
        """Test computing seq stats on taxonomy with only one level."""
        exp = {'1': [3, 4, 'AGGT'], '3': [3, 4, 'AGGC'], '2': [1, 4, 'AGGA']}
        obs = compute_sequence_stats(self.fasta1, self.tax_map3, ['Z'])
        self.assertEqual(obs, exp)

        exp = {'1': [2, 4, 'AGGT'], '3': [4, 4, 'AGGC'], '2': [0, 4, 'AGGA']}
        obs = compute_sequence_stats(self.fasta1, self.tax_map3, ['A'])
        self.assertEqual(obs, exp)

    def test_compute_sequence_stats_multiple_unknown_keywords(self):
        """Test computing seq stats on taxonomy with multiple unknown kws."""
        exp = {'1': [3, 4, 'AGGT'], '3': [2, 4, 'AGGC'], '2': [3, 4, 'AGGA']}
        obs = compute_sequence_stats(self.fasta1, self.tax_map2, ['Z', 'F'])
        self.assertEqual(obs, exp)

    def test_compute_sequence_stats_no_unknowns(self):
        """Test computing seq stats on taxonomy with no unknown kws."""
        exp = {'1': [3, 4, 'AGGT'], '2': [3, 4, 'AGGA'], '3': [5, 4, 'AGGC']}
        obs = compute_sequence_stats(self.fasta1, self.tax_map2)
        self.assertEqual(obs, exp)

        obs = compute_sequence_stats(self.fasta1, self.tax_map2, [])
        self.assertEqual(obs, exp)

    def test_compute_sequence_stats_all_unknown(self):
        """Test computing seq stats on taxonomy with all unknown tax levels."""
        exp = {'1': [0, 4, 'AGGT'], '2': [0, 4, 'AGGA'], '3': [0, 4, 'AGGC']}
        obs = compute_sequence_stats(self.fasta1, self.tax_map2, ['A', 'B',
                'C', 'D', 'F', 'Z'])
        self.assertEqual(obs, exp)

    def test_compute_sequence_stats_invalid_tax_map(self):
        """Test computing seq stats on invalid taxonomy maps."""
        self.assertRaises(ValueError, compute_sequence_stats, self.fasta1,
                          self.tax_map_invalid1, ['Z', 'F'])
        self.assertRaises(ValueError, compute_sequence_stats, self.fasta1,
                          self.tax_map_invalid2, ['Z', 'F'])


if __name__ == "__main__":
    main()
