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

"""Test suite for the summarize_taxonomic_agreement.py module."""

from cogent.util.unit_test import TestCase, main
from nested_reference_otus.summarize_taxonomic_agreement import (
        generate_taxonomic_agreement_summary, parse_taxonomic_information,
        summarize_taxonomic_agreement)

class SummarizeTaxonomicAgreementTests(TestCase):
    """Tests for the summarize_taxonomic_agreement.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.tax_map1 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tG1\tA;B;C;\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;Z;T\tfoo\n"]

        self.otu_map1 = ["A\t1\t2\t3\n"]

        # OTU with only a reference sequence.
        self.otu_map2 = ["A\t1\t2\n", "B\t3\n"]

        # Taxonomy string with no semicolons (only one level).
        self.tax_map2 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tG1\tA\tfoo\n",
                "2\tG2\tB\tfoo\n",
                "3\tG3\tF\tfoo\n"]

        # Taxonomy string with double semicolons and whitespace taxonomy level.
        self.tax_map3 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tG1\tA;B;;C\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;B;;Z;   \tfoo\n"]

        # Missing column in one of the rows.
        self.tax_map_invalid1 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tA;B;C\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;B;Z\tfoo\n"]

        # Missing header.
        self.tax_map_invalid2 = [
                "1\tG1\tA;B;C\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;B;F\tfoo\n"]

        # Wrong number of taxonomic levels.
        self.tax_map_invalid3 = [
                "ID Number\tGenBank Number\tNew Taxon String\tSource\n",
                "1\tG1\tA;B;C\tfoo\n",
                "2\tG2\tA;B;D\tfoo\n",
                "3\tG3\tA;B;Z;Z;F\tfoo\n"]

    def test_parse_taxonomic_information_standard(self):
        """Test parsing a typical taxonomy map."""
        exp = {'1': ['A', 'B', 'C'], '2': ['A', 'B', 'D'],
               '3': ['A', 'Z', 'T']}
        obs = parse_taxonomic_information(self.tax_map1, 3)
        self.assertEqual(obs, exp)

    def test_parse_taxonomic_information_single_level(self):
        """Test parsing a taxonomy map with single-level taxonomy."""
        exp = {'1': ['A'], '2': ['B'], '3': ['F']}
        obs = parse_taxonomic_information(self.tax_map2, 1)
        self.assertEqual(obs, exp)

    def test_parse_taxonomic_information_empty_levels(self):
        """Test parsing a taxonomy map with empty taxonomic levels."""
        exp = {'1': ['A', 'B', 'C'], '2': ['A', 'B', 'D'],
               '3': ['A', 'B', 'Z']}
        obs = parse_taxonomic_information(self.tax_map3, 3)
        self.assertEqual(obs, exp)

    def test_parse_taxonomic_information_invalid_tax_map(self):
        """Test parsing invalid taxonomy maps."""
        self.assertRaises(ValueError, parse_taxonomic_information,
                          self.tax_map_invalid1, 3)
        self.assertRaises(ValueError, parse_taxonomic_information,
                          self.tax_map_invalid2, 3)
        self.assertRaises(ValueError, parse_taxonomic_information,
                          self.tax_map_invalid3, 3)

    def test_generate_taxonomic_agreement_summary_standard(self):
        """Test computing the taxonomic agreement for seqs in OTUs."""
        exp = {'A': [[100.0, 50.0, 0.0], [['A'], ['B', 'Z'], ['C', 'D', 'T']]]}
        obs = generate_taxonomic_agreement_summary(self.otu_map1,
                                                   self.tax_map1, 3)
        self.assertEqual(obs, exp)

    def test_generate_taxonomic_agreement_summary_ref_only(self):
        """Test computing the taxonomic agreement for OTU with only a ref."""
        exp = {'A': [[100.0, 100.0, 0.0], [['A'], ['B'], ['C', 'D']]],
               'B': [[0.0, 0.0, 0.0], [['A'], ['Z'], ['T']]]}
        obs = generate_taxonomic_agreement_summary(self.otu_map2,
                                                   self.tax_map1, 3)
        self.assertEqual(obs, exp)

    def test_summarize_taxonomic_agreement_standard(self):
        """Test writing out the taxonomic agreement for seqs in OTUs."""
        exp = ['A\t100.00\t50.00\t0.00\tA\tB,Z\tC,D,T\n']
        obs = summarize_taxonomic_agreement(self.otu_map1, self.tax_map1, 3)
        self.assertEqual(obs, exp)

    def test_summarize_taxonomic_agreement_multiple_otus(self):
        """Test writing out the taxonomic agreement for multiple OTUs."""
        exp = ['A\t100.00\t100.00\t0.00\tA\tB\tC,D\n',
               'B\t0.00\t0.00\t0.00\tA\tZ\tT\n']
        obs = summarize_taxonomic_agreement(self.otu_map2, self.tax_map1, 3)
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
