#!/usr/bin/env python
# File created on 04 Feb 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
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
from cogent.parse.greengenes import SpecificGreengenesParser
from cogent.app.util import get_tmp_filename
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.workflow import (print_commands,
                            print_to_stdout,
                            no_status_updates,
                            WorkflowLogger,
                            WorkflowError,
                            generate_log_fp)
from qiime.util import (create_dir)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_tree_fp',help='the input tree file (the tips in the tree define which sequences are included from the database)'),
 make_option('-g','--gg_database_fp',help='the greengenes fields fp (gzipped)'),
 make_option('-o','--output_dir',help='the output dir'),
 make_option('-r','--run_id',help='the run id (used in naming some files)'),
 make_option('-s','--similarity_thresholds',help='the similarity thresholds'),
]
script_info['optional_options'] = [
 #make_option('-p','--otu_prefix',help='the otu prefix',default=''),
 make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
 make_option('-w','--print_only',action='store_true',\
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]',default=False),\
]
script_info['version'] = __version__


def call_commands(commands,
                  status_update_callback,
                  logger):
    """Run list of commands, one after another """
    logger.write("Executing commands.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('%s\n%s' % e)
            logger.write('# %s command \n%s\n\n' % e)
            proc = Popen(e[1],shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=PIPE)
            # communicate pulls all stdout/stderr from the PIPEs to 
            # avoid blocking -- don't remove this line!
            stdout, stderr = proc.communicate()
            return_value = proc.returncode
            if return_value != 0:
                msg = "\n\n*** ERROR RAISED DURING STEP: %s\n" % e[0] +\
                 "Command run was:\n %s\n" % e[1] +\
                 "Command returned exit status: %d\n" % return_value +\
                 "Stdout:\n%s\nStderr\n%s\n" % (stdout,stderr)
                logger.write(msg)
                logger.close()
                raise WorkflowError, msg


## Begin task-specific workflow functions
def run_pick_nested_greengenes_otus(input_tree_fp,
                              gg_database_fp,
                              output_dir,
                              run_id,
                              similarity_thresholds,
                              otu_prefix,
                              command_handler,
                              status_update_callback=print_to_stdout):


    # Prepare some variables for the later steps
    create_dir(output_dir)
    otu_dir = join(output_dir,'otus')
    create_dir(otu_dir)
    rep_set_dir = join(output_dir,'rep_set')
    create_dir(rep_set_dir)
    tax_dir = join(output_dir,'taxonomies')
    create_dir(tax_dir)
    tree_dir = join(output_dir,'trees')
    create_dir(tree_dir)
    commands = []
    files_to_remove = []
    
    logger = WorkflowLogger(generate_log_fp(output_dir))
    similarity_thresholds.sort()
    similarity_thresholds.reverse()
    full_tree = LoadTree(input_tree_fp)
    
    gg_data = SpecificGreengenesParser(gzip.open(gg_database_fp),
     fields=['prokMSA_id','decision','aligned_seq','non_ACGT_percent'],
     ids=[tip.Name for tip in full_tree.tips()])
    # create some temp files to store the gg database information
    # that we're interested in
    aligned_input_fp =   get_tmp_filename(prefix='NestedGG',
                                          suffix='.fasta')
    aligned_input_f = open(aligned_input_fp,'w')
    unaligned_input_fp = get_tmp_filename(prefix='NestedGG',
                                          suffix='.fasta')
    unaligned_input_f = open(unaligned_input_fp,'w')
    isolate_fp = get_tmp_filename(prefix='NestedGG',
                                  suffix='.txt')
    isolate_f = open(isolate_fp,'w')
    isolate_f.write('#prokMSA_id\tdecision\tnon_ACGT_percent\taligned_seq\n')
    
    entry_count = 0
    for prokmsa,decision,aligned_seq,non_ACGT_percent in gg_data:
        entry_count += 1
        isolate_f.write('%s\t%s\t%s\tNA\n' % (prokmsa,decision,non_ACGT_percent))
        clean_aligned_seq = aligned_seq.replace('.','-').replace('U','T')
        aligned_input_f.write('>%s\n%s\n' %\
         (prokmsa,clean_aligned_seq))
        unaligned_input_f.write('>%s\n%s\n' %\
         (prokmsa,clean_aligned_seq.replace('-','')))
    logger.write('GG Entries included in analysis: %d\n' % entry_count)
    if not entry_count:
        logger.write("No seqs retained from Greengenes.")
        logger.close()
        raise WorkflowError, "No seqs retained from Greengenes."
    isolate_f.close()
    aligned_input_f.close()
    unaligned_input_f.close()
    
    current_inseqs = unaligned_input_fp
    current_aligned_inseqs = aligned_input_fp
    previous_otu_map = None
    for similarity_threshold in similarity_thresholds:
        # prep presort command
        sorted_fasta_fp = get_tmp_filename(prefix='NestedGG',
                                           suffix='.fasta')
        sorted_fasta_basename = splitext(split(sorted_fasta_fp)[1])[0]
        files_to_remove.append(sorted_fasta_fp)
        presort_cmd = 'presort_gg_fasta.py -f %s -i %s -o %s' % (
                      current_inseqs,
                      isolate_fp,
                      sorted_fasta_fp)
        commands.append([('Presort (%d)' % similarity_threshold, 
                           presort_cmd)])
        
        # pick otus command
        otu_fp = '%s/gg_%d_otu_map.txt' % (otu_dir,similarity_threshold)
        clusters_fp = '%s/gg_%d_clusters.uc' % (otu_dir,similarity_threshold)
        temp_otu_fp = '%s/%s_otus.txt' % (otu_dir, sorted_fasta_basename)
        temp_log_fp = '%s/%s_otus.log' % (otu_dir, sorted_fasta_basename)
        temp_clusters_fp = '%s/%s_clusters.uc' % (otu_dir, sorted_fasta_basename)
        pick_otus_cmd = \
         'pick_otus.py -m uclust -DB -i %s -s %1.2f -o %s' % (
           sorted_fasta_fp,
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
        temp_rep_set_fp = get_tmp_filename(prefix='NestedGG',
                                           suffix='.fasta')
        files_to_remove.append(temp_rep_set_fp)
        pick_rep_set_cmd = \
         'pick_rep_set.py -m first -i %s -o %s -f %s' % (
          otu_fp, 
          temp_rep_set_fp,
          sorted_fasta_fp)
        commands.append([('Pick Rep Set (%d)' % similarity_threshold,
                           pick_rep_set_cmd)])
        
        # Call the command handler on the list of commands
        command_handler(commands, status_update_callback, logger)
        commands = []
        
        # Rename sequences
        rep_set_fp = '%s/gg_%d_otus_%s.fasta' % (
          rep_set_dir,
          similarity_threshold,
          run_id)
        rep_set_f = open(rep_set_fp,'w')
        included_seq_ids = []
        for seq_id, seq in MinimalFastaParser(open(temp_rep_set_fp)):
            seq_id_fields = seq_id.split()
            included_seq_ids.append(seq_id_fields[1])
            rep_set_f.write('>%s\n%s\n' % (' '.join(seq_id_fields[1:]),
                                             seq))
        rep_set_f.close()
        
        # Filter alignment to correspond to rep set
        aligned_rep_set_fp = '%s/gg_%d_otus_%s_aligned.fasta' % (
          rep_set_dir,
          similarity_threshold,
          run_id)
        filter_fasta_cmd = 'filter_fasta.py -f %s -o %s -a %s ' % (
          current_aligned_inseqs,
          aligned_rep_set_fp,
          rep_set_fp)
        commands.append([('Filter aligned seqs (%d)' % similarity_threshold,
                          filter_fasta_cmd)])
                          
        ## Can't currently expand OTU maps as the OTU ids in the OTU
        ## maps are the arbitrary ascending integers.
        # Merge OTU maps so OTUs expand to contain all relvant GG 
        # prokMSA ids (rather than just the previous OTU id)
        # if not previous_otu_map:
        #     previous_otu_map = otu_fp
        # else:
        #     merged_otu_map = \
        #      '%s/gg_%d_merged_otu_map.txt' % (otu_dir,similarity_threshold)
        #     merge_otu_map_cmd = 'merge_otu_maps.py -i %s,%s -o %s' % (
        #      previous_otu_map,
        #      otu_fp,
        #      merged_otu_map)
        #     commands.append([('Merge OTU map (%d)' % similarity_threshold,
        #                        merge_otu_map_cmd)])           
        #     previous_otu_map = merged_otu_map
 
        # Call the command handler on the list of commands
        command_handler(commands, status_update_callback, logger)
        commands = []
        
        # Write tree containing only tips that are in the current 
        # otu set
        tree_fp = '%s/gg_%d_otus_%s.tre' % (
          tree_dir,
          similarity_threshold,
          run_id)
        full_tree.getSubTree(included_seq_ids).writeToFile(tree_fp)
        
        # prep for the next iteration
        remove_files(files_to_remove)
        files_to_remove = []
        current_inseqs = rep_set_fp
        current_aligned_inseqs = aligned_rep_set_fp
        
    logger.close()
    remove_files([unaligned_input_fp,
                  aligned_input_fp,
                  isolate_fp])

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose
    
    input_tree_fp = opts.input_tree_fp
    output_dir = opts.output_dir
    gg_database_fp = opts.gg_database_fp
    run_id = opts.run_id
    similarity_thresholds = map(int,opts.similarity_thresholds.split(','))
    
    otu_prefix = '' # this option is not currently supported opts.otu_prefix
    verbose = opts.verbose
    print_only = opts.print_only
    
    try:
        makedirs(output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            print "Output directory already exists. Please choose "+\
             "a different directory, or force overwrite with -f."
            exit(1)
        
    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands
    
    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    run_pick_nested_greengenes_otus(
     input_tree_fp=input_tree_fp,
     gg_database_fp=gg_database_fp,
     output_dir=output_dir,
     run_id=run_id,
     similarity_thresholds=similarity_thresholds,
     otu_prefix=otu_prefix,
     command_handler=command_handler,
     status_update_callback=status_update_callback)


if __name__ == "__main__":
    main()