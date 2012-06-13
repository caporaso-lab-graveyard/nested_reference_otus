#!/usr/bin/env python
# File created on 04 Feb 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

import gzip
from os import makedirs
from subprocess import Popen, PIPE, STDOUT
from optparse import make_option
from cogent import LoadTree
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.workflow import (print_commands,
                            no_status_updates,
                            WorkflowError,
                            call_commands_serially)

from nested_reference_otus.nested_reference_workflow import (get_second_field,
        rename_rep_seqs, pick_nested_reference_otus)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_fasta_fp',help='the full fasta file to generate reference OTUs from'),
 make_option('-o','--output_dir',help='the output dir'),
 make_option('-r','--run_id',help='the run id (used in naming some files)'),
 make_option('-s','--similarity_thresholds',help='the similarity thresholds'),
]
script_info['optional_options'] = [
 make_option('-w','--print_only',action='store_true',\
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]',default=False),\
 make_option('-t','--input_tree_fp',help='the full tree to filter to otu trees'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    verbose = opts.verbose
    
    input_fasta_fp = opts.input_fasta_fp
    input_tree_fp = opts.input_tree_fp
    output_dir = opts.output_dir
    run_id = opts.run_id
    similarity_thresholds = map(int,opts.similarity_thresholds.split(','))
    
    verbose = opts.verbose
    print_only = opts.print_only
    
    try:
        makedirs(output_dir)
    except OSError:
        print "Output directory already exists. Please choose "+\
         "a different directory, or force overwrite with -f."
        exit(1)
        
    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
    
    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    pick_nested_reference_otus(
     input_fasta_fp=input_fasta_fp,
     input_tree_fp=input_tree_fp,
     output_dir=output_dir,
     run_id=run_id,
     similarity_thresholds=similarity_thresholds,
     command_handler=command_handler,
     status_update_callback=status_update_callback)


if __name__ == "__main__":
    main()
