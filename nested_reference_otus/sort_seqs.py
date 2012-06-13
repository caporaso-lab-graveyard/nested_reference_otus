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

"""Contains functions used in the sort_seqs.py script."""

from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import fields_to_dict

def compute_sequence_stats(fasta_lines, tax_map_lines, unknown_keywords=None):
    seq_stats = {}

    if tax_map_lines[0] != \
            "ID Number\tGenBank Number\tNew Taxon String\tSource\n":
        raise ValueError("The taxonomy map file appears to be invalid "
                         "because it is either missing the header or has a "
                         "corrupt header.")

    # Record the taxonomy depths for each sequence.
    for seq_id, seq_info in fields_to_dict(tax_map_lines[1:]).items():
        if len(seq_info) != 3:
            raise ValueError("The taxonomy map file appears to be invalid "
                             "because it does not have exactly 4 columns.")
        taxonomy = seq_info[1].split(';')

        # Remove any 'unknown' taxonomy levels before computing the known
        # taxonomy depth.
        if unknown_keywords:
            for unknown_keyword in unknown_keywords:
                while unknown_keyword in taxonomy:
                    taxonomy.remove(unknown_keyword)
        seq_stats[seq_id] = [len(taxonomy)]

    # Record the sequence and sequence for each sequence.
    for seq_id, seq in MinimalFastaParser(fasta_lines):
        if seq_id in seq_stats:
            seq_stats[seq_id].extend([len(seq), seq])
        else:
            print ("Found sequence id %s in FASTA file that wasn't in the "
                   "taxonomy mapping file\n" % seq_id)
            # For now, assign a taxonomic depth of 0 because we don't have any
            # taxonomic information for the sequence.
            seq_stats[seq_id] = [0, len(seq), seq]
    return seq_stats
