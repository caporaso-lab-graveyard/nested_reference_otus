#!/usr/bin/env python
from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)
from nested_reference_otus.unnest import parse_otu_map, make_nodes, join_nodes

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Unnest OTUs into a tree"""
script_info['script_description'] = """
Unnest the OTUs to reconstruct the a as to allow a query of all tips that 
are associated with a particular OTU.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Unnests the OTUs",
"Take the nested OTU maps and unroll them into a tree",
"%prog -i gg_99_otu_map.txt,gg_97_otu_map.txt,gg_94_otu_map.txt -o unnested.ntree"))
script_info['output_description']= """
A newick string representing the relationships between the OTU clusters
"""

script_info['required_options'] = [
    make_option('-i','--input_otu_maps',
        help="The input OTU maps. This should be a comma seperated list of "
        "each OTU maps as produced by uclust"),
    options_lookup['output_fp']
]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # expects maps to be in assembly order, ie:
    # gg_99_otu_map.txt,gg_97_otu_map.txt,gg_94_otu_map.txt,...
    otus_in_order = opts.input_otu_maps.split(',')

    last_level = 100 
    nodes = []
    for otus in otus_in_order:
        tmp1,level,tmp2,tmp3 = otus.split('_')

        level = int(level)
        length = float(last_level - level)
        otu_map = parse_otu_map(open(otus))
    
        nodes.append(make_nodes(otu_map, length, level))
    
        last_level = level

    tree = join_nodes(nodes)
    f = open(opts.output_fp,'w')
    f.write(tree.getNewick(with_distances=True))
    f.close()

if __name__ == "__main__":
    main()
