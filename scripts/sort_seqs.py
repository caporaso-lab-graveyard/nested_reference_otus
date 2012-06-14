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

from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)
from nested_reference_otus.sort_seqs import (compute_sequence_stats,
                                             sort_seqs_by_taxonomic_depth)

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

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    seq_stats = compute_sequence_stats(
            open(opts.input_fasta_fp, 'U').readlines(),
            open(opts.input_taxonomy_map, 'U').readlines(),
            ['Incertae_sedis', 'unidentified'])
    seq_stats_sorted = sort_seqs_by_taxonomic_depth(seq_stats)

    # Write out our sorted sequences.
    out_fasta_f = open(opts.output_fp, 'w')
    for seq in seq_stats_sorted:
        out_fasta_f.write('>' + seq[0] + '\n' + seq[3] + '\n')
    out_fasta_f.close()


if __name__ == "__main__":
    main()
