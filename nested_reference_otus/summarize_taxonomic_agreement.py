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

"""Contains functions used in the summarize_taxonomic_agreement.py script."""

from qiime.parse import fields_to_dict

def summarize_taxonomic_agreement(otu_map_lines, tax_map_lines,
                                  taxonomic_levels=8):
    taxonomic_agreement = generate_taxonomic_agreement_summary(otu_map_lines,
            tax_map_lines, taxonomic_levels)
    # Build up our results file in the order of the original OTU map.
    results = []
    for line in otu_map_lines:
        otu_id = line.split('\t')[0]
        agreement_info = taxonomic_agreement[otu_id]

        result_str = '%s' % otu_id
        for level_agreement in agreement_info[0]:
            result_str += '\t%.2f' % level_agreement
        for levels_encountered in agreement_info[1]:
            result_str += '\t' + ','.join(levels_encountered)
        results.append(result_str + '\n')
    return results

def generate_taxonomic_agreement_summary(otu_map_lines, tax_map_lines,
                                         taxonomic_levels=8):
    tax_map = parse_taxonomic_information(tax_map_lines, taxonomic_levels)
    otu_map = fields_to_dict(otu_map_lines)

    taxonomic_agreement = {}
    for otu_id, seq_ids in otu_map.items():
        taxonomic_agreement[otu_id] = [[], []]
        ref_seq_id = seq_ids[0]
        ref_seq_tax = tax_map[ref_seq_id]

        # Calculate percent agreement for each taxonomic level.
        for level_idx, ref_level in enumerate(ref_seq_tax):
            agreement = 0
            encountered_levels = [ref_level]
            for seq_id in seq_ids[1:]:
                seq_level = tax_map[seq_id][level_idx]
                if ref_level == seq_level:
                    agreement += 1
                if seq_level not in encountered_levels:
                    encountered_levels.append(seq_level)
            if len(seq_ids[1:]) > 0:
                taxonomic_agreement[otu_id][0].append(
                        (agreement / len(seq_ids[1:])) * 100)
            else:
                # Prevents divide-by-zero error.
                taxonomic_agreement[otu_id][0].append(0)
            taxonomic_agreement[otu_id][1].append(encountered_levels)
    return taxonomic_agreement

def parse_taxonomic_information(tax_map_lines, taxonomic_levels=8):
    tax_info = {}

    if tax_map_lines[0] != \
            "ID Number\tGenBank Number\tNew Taxon String\tSource\n":
        raise ValueError("The taxonomy map file appears to be invalid "
                         "because it is either missing the header or has a "
                         "corrupt header.")
    for seq_id, seq_info in fields_to_dict(tax_map_lines[1:]).items():
        if len(seq_info) != 3:
            raise ValueError("The taxonomy map file appears to be invalid "
                             "because it does not have exactly 4 columns.")
        # Split at each level and remove any empty levels or levels that
        # contain only whitespace.
        taxonomy = [level for level in seq_info[1].split(';') \
                    if level.strip() != '']
        if len(taxonomy) != taxonomic_levels:
            raise ValueError("Encountered invalid taxonomy '%s'. Valid "
                    "taxonomy strings must have %d levels separated by "
                    "semicolons." % (seq_info[1], taxonomic_levels))
        tax_info[seq_id] = taxonomy
    return tax_info
