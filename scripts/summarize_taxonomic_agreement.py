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
from nested_reference_otus.summarize_taxonomic_agreement import (
        summarize_taxonomic_agreement)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Summarizes the agreement of sequences with their reference taxonomy"
script_info['script_description'] = """
This script summarizes the agreement of each sequence's taxonomy with the
reference sequence taxonomy in the OTU that the sequence resides in. The
taxonomic agreement is expressed as a percentage at each taxonomic level. The
script also reports all values that are encountered at a particular taxonomic
level, with the reference's value listed first.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Summarize taxonomic agreement",
"Summarizes the percentage of sequences in each OTU that have the same "
"taxonomic level as the reference sequence", "%prog -i 99_otu_map.txt -t "
"taxonomy_map.txt -o taxonomic_agreement_summary.txt"))
script_info['output_description']= """
The script creates a single tab-separated file containing the taxonomic
agreement summary.
"""

script_info['required_options'] = [
    options_lookup['otu_map_as_primary_input'],
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

    results = summarize_taxonomic_agreement(
            open(opts.otu_map_fp, 'U').readlines(),
            open(opts.input_taxonomy_map, 'U').readlines())

    out_f = open(opts.output_fp, 'w')
    out_f.write('OTU_ID\tDomain\tKingdom\tPhylum\tClass\tOrder\tFamily\t'
                'Genus\tSpecies\tDomain\tKingdom\tPhylum\tClass\tOrder\t'
                'Family\tGenus\tSpecies\n')
    for line in results:
        out_f.write(line)
    out_f.close()


if __name__ == "__main__":
    main()
