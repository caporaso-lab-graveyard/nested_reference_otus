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

from StringIO import StringIO
import sys
from cogent.util.unit_test import TestCase, main
from nested_reference_otus.sort_seqs import (compute_sequence_stats,
                                             sort_seqs_by_taxonomic_depth)

class SortSeqsTests(TestCase):
    """Tests for the sort_seqs.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        # Equal read lengths.
        self.fasta1 = [">1", "AGGT", ">2", "AGGA", ">3", "AGGC"]

        # Unequal read lengths.
        self.fasta2 = [">1", "AGGTAC", ">2", "AG", ">3", "AGGCAAA"]

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

        # Sorting keys are all the same.
        self.seq_stats1 = {'1': [3, 4, 'AGGT'], '3': [3, 4, 'AGGC'],
                           '2': [3, 4, 'AGGA']}

        # All different keys.
        self.seq_stats2 = {'1': [3, 4, 'AGGT'], '3': [5, 1, 'A'],
                           '2': [6, 7, 'AGGAGTC']}

        # Single sequence.
        self.seq_stats3 = {'4': [3, 4, 'AGGT']}

        # Some taxonomic depths the same.
        self.seq_stats4 = {'1': [3, 4, 'AGGT'], '3': [3, 5, 'AGGTA'],
                           '2': [6, 7, 'AGGAGTC']}

        # All taxonomic depths the same.
        self.seq_stats5 = {'1': [3, 4, 'AGGT'], '3': [3, 5, 'AGGTA'],
                           '2': [3, 7, 'AGGAGTC']}

        # Missing sequence in FASTA, but found in taxonomy mapping file.
        self.seq_stats6 = {'42': [8], '3': [3, 5, 'AGGTA'],
                           '2': [3, 7, 'AGGAGTC']}


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

    def test_compute_sequence_stats_unequal_read_lengths(self):
        """Test computing seq stats on unequal read lengths."""
        exp = {'1': [3, 6, 'AGGTAC'], '3': [2, 7, 'AGGCAAA'],
               '2': [3, 2, 'AG']}
        obs = compute_sequence_stats(self.fasta2, self.tax_map1, ['Z'])
        self.assertEqual(obs, exp)

    def test_compute_sequence_stats_missing_taxonomy(self):
        """Test computing seq stats on a sequence with no taxonomy mapping."""
        # Save stdout and replace it with something that will capture the print
        # statement. Note: this code was taken from here:
        # http://stackoverflow.com/questions/4219717/how-to-assert-output-
        #     with-nosetest-unittest-in-python/4220278#4220278
        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out

            self.fasta1.append(">4")
            self.fasta1.append("AACCGGTT")
            exp = {'1': [3, 4, 'AGGT'], '3': [2, 4, 'AGGC'], '2': [3, 4, 'AGGA'],
                   '4': [0, 8, 'AACCGGTT']}
            obs = compute_sequence_stats(self.fasta1, self.tax_map1, ['Z'])
            self.assertEqual(obs, exp)

            output = out.getvalue().strip()
            self.assertEqual(output, "Found sequence id '4' in the FASTA file "
                             "that wasn't in the taxonomy mapping file")
        finally:
            sys.stdout = saved_stdout

    def test_compute_sequence_stats_invalid_tax_map(self):
        """Test computing seq stats on invalid taxonomy maps."""
        self.assertRaises(ValueError, compute_sequence_stats, self.fasta1,
                          self.tax_map_invalid1, ['Z', 'F'])
        self.assertRaises(ValueError, compute_sequence_stats, self.fasta1,
                          self.tax_map_invalid2, ['Z', 'F'])

    def test_sort_seqs_by_taxonomic_depth_equal_keys(self):
        """Test sorting sequences when keys are all equal."""
        exp = [['1', 3, 4, 'AGGT'], ['3', 3, 4, 'AGGC'], ['2', 3, 4, 'AGGA']]
        obs = sort_seqs_by_taxonomic_depth(self.seq_stats1)
        # Since we were given a dictionary of sequence stats, the order will
        # stay the same as was in the dictionary since we are doing a stable
        # sort. Thus, we can't simply test for equality because the key order
        # in the dictionary is not guaranteed.
        self.assertEqual(obs[0][1:3], exp[0][1:3])
        self.assertEqual(obs[1][1:3], exp[1][1:3])
        self.assertEqual(obs[2][1:3], exp[2][1:3])

    def test_sort_seqs_by_taxonomic_depth_unequal_keys(self):
        """Test sorting sequences when keys are all unequal."""
        exp = [['2', 6, 7, 'AGGAGTC'], ['3', 5, 1, 'A'], ['1', 3, 4, 'AGGT']]
        obs = sort_seqs_by_taxonomic_depth(self.seq_stats2)
        self.assertEqual(obs, exp)

    def test_sort_seqs_by_taxonomic_depth_single_seq(self):
        """Test sorting a single sequence."""
        exp = [['4', 3, 4, 'AGGT']]
        obs = sort_seqs_by_taxonomic_depth(self.seq_stats3)
        self.assertEqual(obs, exp)

    def test_sort_seqs_by_taxonomic_depth_same_tax_depth(self):
        """Test sorting with some taxonomic depths being the same."""
        exp = [['2', 6, 7, 'AGGAGTC'], ['3', 3, 5, 'AGGTA'],
               ['1', 3, 4, 'AGGT']]
        obs = sort_seqs_by_taxonomic_depth(self.seq_stats4)
        self.assertEqual(obs, exp)

    def test_sort_seqs_by_taxonomic_depth_all_same_tax_depth(self):
        """Test sorting with all taxonomic depths being the same."""
        exp = [['2', 3, 7, 'AGGAGTC'], ['3', 3, 5, 'AGGTA'],
               ['1', 3, 4, 'AGGT']]
        obs = sort_seqs_by_taxonomic_depth(self.seq_stats5)
        self.assertEqual(obs, exp)

    def test_sort_seqs_by_taxonomic_depth_mising_seq_in_fasta(self):
        """Test sorting with a sequence that was missing from FASTA file."""
        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out

            exp = [['2', 3, 7, 'AGGAGTC'], ['3', 3, 5, 'AGGTA']]
            obs = sort_seqs_by_taxonomic_depth(self.seq_stats6)
            self.assertEqual(obs, exp)

            output = out.getvalue().strip()
            self.assertEqual(output, "Found sequence id '42' in the taxonomy "
                    "mapping file that wasn't in the FASTA file")
        finally:
            sys.stdout = saved_stdout


if __name__ == "__main__":
    main()
