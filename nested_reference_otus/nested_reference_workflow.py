#!/usr/bin/env python
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


"""Contains functions used in the nested_reference_workflow.py script."""

from os.path import join, split, splitext
from cogent.app.util import get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import remove_files
from qiime.util import create_dir
from qiime.workflow.util import generate_log_fp, print_to_stdout, WorkflowLogger

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
         'pick_otus.py -m uclust -DBz -i %s -s %1.2f -o %s' % (
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
