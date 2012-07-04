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

"""Contains functions used in the summarize_taxonomic_agreement.py script."""

from qiime.parse import fields_to_dict

def summarize_taxonomic_agreement(otu_map_lines, tax_map_lines,
                                  taxonomic_levels=8):
    """Computes a summary of taxonomic agreement between ref and its seqs.

    Returns a list of lines suitable for writing to an output file. Each line
    is for a single OTU, and the OTUs are ordered the same as the input OTU
    map. The columns are separated by tabs.
    
    The first column is the OTU ID, followed by the OTU size (number of seqs,
    including the ref), followed by a comma-separated list of sequence
    identifiers for all sequences in the OTU (including the reference, which is
    always listed first). Next, there is a column for each taxonomic level
    found in the input taxonomy mapping file. Each column contains the percent
    agreement at that particular level. After that, additional columns are
    provided for each taxonomic level that contain a comma-separated list of
    taxonomic values that were encountered at that level, with the reference
    taxonomy listed first.

    Arguments:
        otu_map_lines - list of lines in the OTU map (the result of calling
            readlines() on the open file handle)
        tax_map_lines - list of lines from the taxonomy mapping file (the
            result of calling readlines() on the open file handle)
        taxonomic_levels - the number of taxonomic levels in the taxonomy
            strings found in the taxonomy mapping file. All taxonomy strings
            must have this number of levels to prevent inconsistent results in
            the summary
    """
    taxonomic_agreement = _generate_taxonomic_agreement_summary(otu_map_lines,
            tax_map_lines, taxonomic_levels)
    # Build up our results file in the order of the original OTU map.
    results = []
    for line in otu_map_lines:
        otu_id = line.split('\t')[0]
        agreement_info = taxonomic_agreement[otu_id]

        result_str = '%s\t%d\t' % (otu_id, agreement_info[0])

        # We put explicit quotes around each seq ID since some of the IDs are
        # numbers and this messes with programs like Excel, where they try to
        # interpret them as a number with commas in it.
        for idx, seq_id in enumerate(agreement_info[1]):
            result_str += "'" + seq_id + "'"
            if idx != len(agreement_info[1]) - 1:
                result_str += ','
        for level_agreement in agreement_info[2]:
            result_str += '\t%.2f%%' % level_agreement
        for levels_encountered in agreement_info[3]:
            result_str += '\t' + ','.join(levels_encountered)
        results.append(result_str + '\n')
    return results

def _generate_taxonomic_agreement_summary(otu_map_lines, tax_map_lines,
                                         taxonomic_levels=8):
    """Computes a summary of taxonomic agreement between ref and its seqs.

    Returns a dictionary with OTU ID as the key. The value is a four-element
    list. The first element is the size of the OTU (i.e. the number of seqs in
    the OTU, including the reference). The second element is a list of sequence
    identifiers for each sequence in the OTU. The reference sequence ID will
    always be listed first, followed by the sequence IDs of the other members
    of the OTU as they appear in the OTU map. The third element is a list
    containing percent agreement at each taxonomic level (a list of floats).
    The fourth element is a list containing all taxonomic values that were
    encountered at each level. The reference taxonomic value will always be
    listed first. The third and fourth elements of the top-level list will
    always be the same length (taxonomic_levels) because they each contain
    information for each taxonomic level.

    Arguments:
        otu_map_lines - list of lines in the OTU map (the result of calling
            readlines() on the open file handle)
        tax_map_lines - list of lines from the taxonomy mapping file (the
            result of calling readlines() on the open file handle)
        taxonomic_levels - the number of taxonomic levels in the taxonomy
            strings found in the taxonomy mapping file. All taxonomy strings
            must have this number of levels to prevent inconsistent results in
            the summary
    """
    tax_map = _parse_taxonomic_information(tax_map_lines, taxonomic_levels)
    otu_map = fields_to_dict(otu_map_lines)

    taxonomic_agreement = {}
    for otu_id, seq_ids in otu_map.items():
        otu_size = len(seq_ids)
        taxonomic_agreement[otu_id] = [otu_size, seq_ids, [], []]

        # The reference sequence is always the first sequence listed in the OTU
        # map.
        ref_seq_id = seq_ids[0]
        ref_seq_tax = tax_map[ref_seq_id]

        # Calculate percent agreement for each taxonomic level. If the OTU only
        # contains a reference sequence, the percent agreement will be 100%.
        # Also keep track of all unique taxonomic values that are encountered
        # for each level (with the reference's taxonomic value listed first).
        for level_idx, ref_level in enumerate(ref_seq_tax):
            agreement_count = 0
            encountered_levels = []
            for seq_id in seq_ids:
                seq_level = tax_map[seq_id][level_idx]
                if ref_level == seq_level:
                    agreement_count += 1
                if seq_level not in encountered_levels:
                    encountered_levels.append(seq_level)
            taxonomic_agreement[otu_id][2].append(
                    (agreement_count / otu_size) * 100)
            taxonomic_agreement[otu_id][3].append(encountered_levels)
    return taxonomic_agreement

def _parse_taxonomic_information(tax_map_lines, taxonomic_levels=8):
    """Parses a taxonomy mapping file to return mapping of seq ID to taxonomy.
    
    Returns a dictionary with sequence ID as the key and a list containing the
    taxonomy at each level. Empty taxonomic levels (i.e. ';;' or levels
    containing only whitespace) are ignored.

    Arguments:
        tax_map_lines - list of lines from the taxonomy mapping file (the
            result of calling readlines() on the open file handle)
        taxonomic_levels - the number of taxonomic levels in the taxonomy
            strings found in the taxonomy mapping file. All taxonomy strings
            must have this number of levels (excluding empty taxonomic levels)
    """
    tax_info = {}

    if tax_map_lines[0] != \
            "ID Number\tGenBank Number\tNew Taxon String\tSource\n":
        raise ValueError("The taxonomy map file appears to be invalid "
                         "because it is either missing the header or has a "
                         "corrupt header.")
    for seq_id, seq_info in fields_to_dict(tax_map_lines[1:]).items():
        if len(seq_info) != 3:
            raise ValueError("The taxonomy map file appears to be invalid "
                             "because it does not have exactly 4 columns.")
        # Split at each level and remove any empty levels or levels that
        # contain only whitespace.
        taxonomy = [level for level in seq_info[1].split(';') \
                    if level.strip() != '']
        if len(taxonomy) != taxonomic_levels:
            raise ValueError("Encountered invalid taxonomy '%s'. Valid "
                    "taxonomy strings must have %d levels separated by "
                    "semicolons." % (seq_info[1], taxonomic_levels))
        tax_info[seq_id] = taxonomy
    return tax_info
