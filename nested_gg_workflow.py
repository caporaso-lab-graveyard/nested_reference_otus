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
 make_option('-i','--input_fasta_fp',help='the full fasta file to generate reference OTUs from'),
 make_option('-o','--output_dir',help='the output dir'),
 make_option('-r','--run_id',help='the run id (used in naming some files)'),
 make_option('-s','--similarity_thresholds',help='the similarity thresholds'),
]
script_info['optional_options'] = [
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
def run_pick_nested_greengenes_otus(input_fasta_fp,
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
    # currently not doing anything with taxonomies and trees
    # tax_dir = join(output_dir,'taxonomies')
    # create_dir(tax_dir)
    # tree_dir = join(output_dir,'trees')
    # create_dir(tree_dir)
    commands = []
    files_to_remove = []
    
    logger = WorkflowLogger(generate_log_fp(output_dir))
    similarity_thresholds.sort()
    similarity_thresholds.reverse()
    
    current_inseqs_fp = input_fasta_fp
    current_inseqs_basename = splitext(split(current_inseqs_fp)[1])[0]
    previous_otu_map = None
    for similarity_threshold in similarity_thresholds:
        
        # pick otus command
        otu_fp = '%s/gg_%d_otu_map.txt' % (otu_dir,similarity_threshold)
        clusters_fp = '%s/gg_%d_clusters.uc' % (otu_dir,similarity_threshold)
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
        temp_rep_set_fp = get_tmp_filename(prefix='NestedGG',
                                           suffix='.fasta')
        files_to_remove.append(temp_rep_set_fp)
        pick_rep_set_cmd = \
         'pick_rep_set.py -m first -i %s -o %s -f %s' % (
          otu_fp, 
          temp_rep_set_fp,
          current_inseqs_fp)
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
        
        
        # prep for the next iteration
        remove_files(files_to_remove)
        files_to_remove = []
        current_inseqs = rep_set_fp
        
    logger.close()

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose
    
    input_fasta_fp = opts.input_fasta_fp
    output_dir = opts.output_dir
    run_id = opts.run_id
    similarity_thresholds = map(int,opts.similarity_thresholds.split(','))
    
    otu_prefix = '' # this option is not currently supported opts.otu_prefix
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
        command_handler = call_commands
    
    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    run_pick_nested_greengenes_otus(
     input_fasta_fp=input_fasta_fp,
     output_dir=output_dir,
     run_id=run_id,
     similarity_thresholds=similarity_thresholds,
     otu_prefix=otu_prefix,
     command_handler=command_handler,
     status_update_callback=status_update_callback)


if __name__ == "__main__":
    main()