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
from os.path import split, splitext, join
from os import makedirs
from subprocess import Popen, PIPE, STDOUT
from optparse import make_option
from cogent import LoadTree
from cogent.util.misc import remove_files
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.util import get_tmp_filename
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.workflow import (print_commands,
                            print_to_stdout,
                            no_status_updates,
                            WorkflowLogger,
                            WorkflowError,
                            generate_log_fp,
                            call_commands_serially)
from qiime.util import (create_dir)

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

def get_second_field(s):
    return s.split()[1]

def rename_rep_seqs(inseqs,rename_f=get_second_field):
    """ """
    for seq_id, seq in MinimalFastaParser(inseqs):
        yield rename_f(seq_id), seq

## Begin task-specific workflow functions
def pick_nested_reference_otus(input_fasta_fp,
                              input_tree_fp,
                              output_dir,
                              run_id,
                              similarity_thresholds,
                              command_handler,
                              status_update_callback=print_to_stdout):


    # Prepare some variables for the later steps
    create_dir(output_dir)
    otu_dir = join(output_dir,'otus')
    create_dir(otu_dir)
    rep_set_dir = join(output_dir,'rep_set')
    create_dir(rep_set_dir)
    # currently not doing anything with taxonomies and trees
    # tax_dir = join(output_dir,'taxonomies')
    # create_dir(tax_dir)
    if input_tree_fp:
        tree_dir = join(output_dir,'trees')
        create_dir(tree_dir)
    commands = []
    files_to_remove = []
    
    logger = WorkflowLogger(generate_log_fp(output_dir))
    similarity_thresholds.sort()
    similarity_thresholds.reverse()
    
    current_inseqs_fp = input_fasta_fp
    current_tree_fp = input_tree_fp
    previous_otu_map = None
    for similarity_threshold in similarity_thresholds:
        current_inseqs_basename = splitext(split(current_inseqs_fp)[1])[0]
        
        # pick otus command
        otu_fp = '%s/%d_otu_map.txt' % (otu_dir,similarity_threshold)
        clusters_fp = '%s/%d_clusters.uc' % (otu_dir,similarity_threshold)
        temp_otu_fp = '%s/%s_otus.txt' % (otu_dir, current_inseqs_basename)
        temp_log_fp = '%s/%s_otus.log' % (otu_dir, current_inseqs_basename)
        temp_clusters_fp = '%s/%s_clusters.uc' % (otu_dir, current_inseqs_basename)
        pick_otus_cmd = \
         'pick_otus.py -m uclust -DB -i %s -s %1.2f -o %s' % (
           current_inseqs_fp,
           similarity_threshold/100,
           otu_dir)
        
        commands.append([('Pick OTUs (%d)' % similarity_threshold,
                          pick_otus_cmd)])
        commands.append([('Rename OTU file (%d)' % similarity_threshold,
                          'mv %s %s' % (temp_otu_fp,otu_fp))])
        commands.append([('Rename uc file (%d)' % similarity_threshold,
                          'mv %s %s' % (temp_clusters_fp,clusters_fp))])
        files_to_remove.append(temp_log_fp)
        
        # rep set picking
        temp_rep_set_fp = get_tmp_filename(prefix='NestedReference',
                                           suffix='.fasta')
        pick_rep_set_cmd = \
         'pick_rep_set.py -m first -i %s -o %s -f %s' % (
          otu_fp, 
          temp_rep_set_fp,
          current_inseqs_fp)
        commands.append([('Pick Rep Set (%d)' % similarity_threshold,
                           pick_rep_set_cmd)])
        command_handler(commands, status_update_callback, logger, close_logger_on_success=False)
        commands = []
        
        # rename representative sequences
        rep_set_fp = '%s/%d_otus_%s.fasta' % (
          rep_set_dir,
          similarity_threshold,
          run_id)
        logger.write('Renaming OTU representative sequences so OTU ids are reference sequence ids.')
        rep_set_f = open(rep_set_fp,'w')
        for e in rename_rep_seqs(open(temp_rep_set_fp,'U')):
            rep_set_f.write('>%s\n%s\n' % e)
        rep_set_f.close()
        files_to_remove.append(temp_rep_set_fp)
        
        # filter the tree, if provided
        if current_tree_fp != None:
            tree_fp = '%s/%d_otus_%s.tre' % (
              tree_dir,
              similarity_threshold,
              run_id)
            tree_cmd = 'filter_tree.py -i %s -f %s -o %s' %\
               (current_tree_fp,rep_set_fp,tree_fp)
            commands.append([('Filter tree (%d)' % similarity_threshold,tree_cmd)])
            command_handler(commands, status_update_callback, logger, close_logger_on_success=False)
            # prep for the next iteration
            current_tree_fp = tree_fp
        
        
        # prep for the next iteration
        remove_files(files_to_remove)
        commands = []
        files_to_remove = []
        current_inseqs_fp = rep_set_fp
        
    logger.close()

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