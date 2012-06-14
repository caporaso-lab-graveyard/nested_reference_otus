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

from operator import itemgetter
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import fields_to_dict

def compute_sequence_stats(fasta_lines, tax_map_lines, unknown_keywords=None):
    """Generates statistics for the input sequences.

    Returns a dictionary with sequence ID as the key, and a list of statistics
    as the value. The statistics (in order of placement in the list) are:
        relevant taxonomic depth (integer)
        sequence read length (integer)
        sequence data (string)
    
    For example, the relevant taxonomic depth for the taxonomy string 'A;B;C'
    would be 3. If unknown_keywords is supplied, any taxonomic level matching
    a value in unknown_keywords wil be ignored. Empty or whitespace-only
    taxonomic levels will also be ignored in the count. For example, if 'Z' is
    an unknown keyword, the relevant taxonomic depth of 'A;B;Z;C' will be 3.

    Arguments:
        fasta_lines - list of lines in FASTA format (the result of calling
            readlines() on the open file handle)
        tax_map_lines - list of lines from the taxonomy mapping file (the
            result of calling readlines() on the open file handle)
        unknown_keywords - a list of strings corresponding to taxonomic level
            strings that should be ignored when computing the relevant
            taxonomic depth
    """
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
        # Split at each level and remove any empty levels or levels that
        # contain only whitespace.
        taxonomy = [level for level in seq_info[1].split(';') \
                    if level.strip() != '']

        # Remove any 'unknown' taxonomy levels before computing the known
        # taxonomy depth.
        if unknown_keywords:
            for unknown_keyword in unknown_keywords:
                while unknown_keyword in taxonomy:
                    taxonomy.remove(unknown_keyword)
        seq_stats[seq_id] = [len(taxonomy)]

    # Record the sequence data and sequence length for each sequence.
    for seq_id, seq in MinimalFastaParser(fasta_lines):
        if seq_id in seq_stats:
            seq_stats[seq_id].extend([len(seq), seq])
        else:
            print ("Found sequence id '%s' in the FASTA file that wasn't in "
                   "the taxonomy mapping file\n" % seq_id)
            # Assign a taxonomic depth of 0 because we don't have any
            # taxonomic information for the sequence.
            seq_stats[seq_id] = [0, len(seq), seq]
    return seq_stats

def sort_seqs_by_taxonomic_depth(seq_stats):
    """Sorts the input sequences by relevant taxonomic depth and read length.

    The sequences are sorted by relevant taxonomic depth (descending) and for
    sequences with the same taxonomic depth, the sequences are sorted by read
    length (descending).

    Returns a list of lists containing the following inner elements:
        sequence ID (string)
        relevant taxonomic depth (integer)
        sequence read length (integer)
        sequence data (string)

    Arguments:
        seq_stats - the output of compute_sequence_stats, which is a dictionary
            mapping sequence ID to sequence statistical information
    """
    # Build a list from our sequence statistics so that we can sort it.
    seq_stats_list = []
    for seq_id, stats in seq_stats.items():
        if len(stats) == 3:
            seq_stats_list.append([seq_id] + stats)
        else:
            print ("Found sequence id '%s' in the taxonomy mapping file that "
                   "wasn't in the FASTA file\n" % seq_id)

    # Sort on sequence length (decreasing) first, then on taxonomic depth (also
    # decreasing). sorted() is guaranteed to give us a stable sort, so we will
    # ultimately end up with the sequences ordered with the most taxonomic
    # information first, with sequences having the same level of taxonomic
    # information ordered by decreasing length.
    return sorted(seq_stats_list, key=itemgetter(1, 2), reverse=True)
