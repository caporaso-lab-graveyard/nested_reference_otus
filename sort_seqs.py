#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

from operator import itemgetter
from cogent.parse.fasta import MinimalFastaParser

from qiime.parse import fields_to_dict
from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Sorts sequences by taxonomic depth and then by length"""
script_info['script_description'] = """
This script sorts a set of input sequences by taxonomic depth (i.e. the amount
of taxonomic information available for a sequence) and by decreasing sequence
read length for sequences that have the same taxonomic depth. This script may
be used to presort sequences before picking OTUs to ensure that OTU centroids
are "quality" sequences.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Sort sequences by taxonomic depth and "
"length",
"Sorts the sequences contained in the input FASTA file by taxonomic depth and "
"then sorts sequences within each taxonomic depth by decreasing length so "
"that longer reads come first.",
"%prog -i input.fasta -t taxonomy_map.txt -o sorted_seqs.fasta"))
script_info['output_description']= """
The script creates a single output FASTA file containing the sorted sequences.
"""

script_info['required_options'] = [
    options_lookup['fasta_as_primary_input'],
    make_option('-t','--input_taxonomy_map',
        help='the input taxonomy map file. This should be a tab-separated '
        'file where each row contains a sequence ID, GenBank number, taxonomy '
        'string, and source string separated by a tab'),
    options_lookup['output_fp']
]
script_info['optional_options'] = []
script_info['version'] = __version__

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

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    seq_stats = compute_sequence_stats(
            open(opts.input_fasta_fp, 'U').readlines(),
            open(opts.input_taxonomy_map, 'U').readlines(),
            ['Incertae_sedis', 'unidentified'])
    seq_stats_sorted = []
    for seq_id, stats in seq_stats.items():
        if len(stats) == 3:
            seq_stats_sorted.append([seq_id] + stats)
        else:
            print ("Found sequence id %s in taxonomy mapping file that wasn't "
                   "in FASTA file\n" % seq_id)

    # Sort on sequence length (decreasing) first, then on taxonomic depth (also
    # decreasing). sorted() is guaranteed to give us a stable sort.
    seq_stats_sorted = sorted(seq_stats_sorted, key=itemgetter(2, 1),
                              reverse=True)
    out_fasta_f = open(opts.output_fp, 'w')
    for seq in seq_stats_sorted:
        out_fasta_f.write('>' + seq[0] + '\n' + seq[3] + '\n')
    out_fasta_f.close()


if __name__ == "__main__":
    main()
