#!/usr/bin/env python
# coding: utf-8

import argparse
from Bio import pairwise2
from glob import glob
from math import sqrt
import numpy as np
import os
import pickle
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import rmsd_atoms
from pyrosetta.rosetta.core.scoring import superimpose_pose
from pyrosetta.rosetta.core.select.residue_selector import     AndResidueSelector, ChainSelector, NotResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import     PerResidueRMSDMetric
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import     OperateOnResidueSubset, RestrictAbsentCanonicalAASRLT
from pyrosetta.rosetta.protocols.grafting import CCDEndsGraftMover
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.util import ChainbreakUtil
from random import choice


# General utility functions

def get_distance(c1, c2):
    """ Returns the distance between two Rosetta XYZ coordinate vectors"""
    return np.linalg.norm(np.array(c1) - np.array(c2))


def find_res_coords(pose, resnum, atom_type='CA'):
    """ For a given pose and residue number, returns the coordinates of a 
    specified atom type """
    residue = pose.residue(resnum)
    res = residue.atom(atom_type)
    return res.xyz()


def list_pose_coords(pose, atom_type='CA'):
    """ For a given pose, list all CA coordinates """
    # Initialize list of CA coordinates
    coords = []

    # Populating list of CA coordinates 
    for res in range(1, pose.total_residue() + 1):
        res_ca = find_res_coords(pose, res, atom_type=atom_type)
        coords.append(res_ca)

    return coords


def check_pose_continuity(pose):
    """
    Scans through all residues in a pose, checking CA and N coordinates, and 
    identifies chain breaks by distances deviating more than 10% from the ideal 
    1.33 A. (This is how the ChainBreakFilter worked.) Returns a bool indicating 
    if loop is continuous (so True means no breaks were found), a list of 
    C-N distances, and a list of pose numbers at which breaks were found. 
    """
    # Get lists of N and C residues
    n_coords = list_pose_coords(pose, atom_type='N')
    c_coords = list_pose_coords(pose, atom_type='C')

    # Check C-N diatances
    continuous = True
    c_n_distances = []
    break_sites = []
    for i in range(len(n_coords) - 1):
        distance = get_distance(c_coords[i], n_coords[i+1])
        c_n_distances.append(distance)

    # Check whether distance indicates a chain break
        if not 0.9 * 1.33 <= distance <= 1.1 * 1.33:
            continuous = False 
            break_sites.append(i)

    return continuous, c_n_distances, break_sites


def selector_to_list(pose, selector):
    """
    Produces a list of residues from a pose identified by a given selector
    """
    from pyrosetta.rosetta.core.select import get_residues_from_subset

    return list(get_residues_from_subset(selector.apply(pose)))


def align_protein_sections(pose_1, selector_1, pose_2, selector_2, mode='CA'):
    """
    Aligns selected regions of two poses, superimposing the second pose onto 
    the first, based on RMSD. Returns the RMSD value. Input selections can
    either be Rosetta selector objects or lists of residue by pose number, 
    which will be made into a selector. By default, uses CA RMSD. Can also 
    use full backbone.
    """
    # Put residue inputs into selector format
    # If a list or range is given, make a selector
    if type(selector_1) in [list, range]:
        select_str_1 = ','.join([str(i) for i in selector_1])
        selector_1 = ResidueIndexSelector(select_str_1)
    if type(selector_2) in [list, range]:
        select_str_2 = ','.join([str(i) for i in selector_2])
        selector_2 = ResidueIndexSelector(select_str_2)

    # Verify mode is acceptable
    assert mode in ['CA', 'BB']
    
    # Set RMSD type based on mode
    if mode == 'CA':
        rmsd_type = rmsd_atoms.rmsd_protein_bb_ca
    if mode == 'BB':
        rmsd_type = rmsd_atoms.rmsd_protein_bb_heavy

    # Set up RMSD metric to align poses
    prmsd = PerResidueRMSDMetric()
    prmsd.set_rmsd_type(rmsd_type)
    prmsd.set_comparison_pose(pose_1)
    prmsd.set_residue_selector_reference(selector_1)
    prmsd.set_residue_selector(selector_2)
    amap = prmsd.create_atom_id_map(pose_2)

    return superimpose_pose(pose_2, pose_1, amap)


def get_b_factor(pose, residue):
    """ 
    Given a pose and a residue number, will return the average b-factor of the 
    backbone atoms (N, CA, C) for the specified residue. Requires residue to  
    be input as a pose number, as opposed to a PDB number. 
    """
    bfactor = pose.pdb_info().bfactor
    atom_index = pose.residue(residue).atom_index

    total_b = 0.0
    for atom in ['N', 'CA', 'C']:
        total_b += bfactor(residue, atom_index(atom))

    # Return average for three atoms
    return total_b / 3


def variable_sliding_window(inp, min_size=0, max_size=0):
    """
    Takes a string or list input and returns a list of frames in that input. The 
    frame size increases with each iteration. Thus, running on an input of 'str' 
    with produce the output ['s', 't', 'r', 'st', 'tr', 'str']. Default behavior 
    will go from a frame size of 1 up to the full size of the input. Optionally, 
    this can be constrained to a set window size range. 
    """
    # Initialize output list
    out_windows_list = []

    # Set initial window size
    if min_size:
        window_size = min_size
    else:
        window_size = 1

    # Set final window size
    if max_size:
        window_max = max_size
    else:
        window_max = len(inp) + 1

    # Outer loop with increasing window size
    while window_size <= window_max:
        frame_start = 0
        frame_end = frame_start + window_size

        # Inner loop sliding the window
        while frame_end <= len(inp):
            # Add frame to output list
            out_windows_list.append(inp[frame_start:frame_end])

            # Increment start of frame and end of frame
            frame_start += 1
            frame_end = frame_start + window_size

        # Increment window size
        window_size += 1

    return out_windows_list


def apply_mask(start_list, mask_list, mode=True):
    """
    Given a list and a masking list, returns a portion of the list as determined 
    by the masking list. By default, returns the members of the list that are 
    identified in the mask as True, but this can be changed to isolate the false 
    members instead.
    """
    # Warning if lengths don't match
    if len(start_list) != len(mask_list):
        print('Warning: List and mask are not the same length.')
        print('         Information may be lost.')

    # Zip list and mask
    phantom = zip(start_list, mask_list)

    if mode:
        return [list_item for list_item,mask_state in phantom if mask_state] 

    else:
        return [list_item for list_item,mask_state in phantom if not mask_state]


# Functions used by protease_info class

class protease_info():
    """
    Data storage structure that encompasses the relationship between two 
    proteins based on their Dali alignments and structural info calculated 
    using Rosetta. The query is assumed to be the reference protein, and the 
    subject a match identified by Dali search. Included are their names, the 
    comparison statistics collected by Dali (Z_SCORE, RMSD, LALI, NRES, PID, 
    and PDB description of the subject), a list of aligned_residue objects 
    corresponding to each residue in the alignment, an identification of 
    catalytic triad residues in the subject, and a list of matched_loop objects 
    reflecting structural comparisons between identified loop regions of the 
    query and subject, which might be exchanged in design applications.
    """
    def __init__(self, query_name, subject_name, dali_file, align_file, 
        query_pose, subject_pose, catalytic_residues, structure_map,  
        auto=True, report=False, verbose=False):

        # Protein names
        self.query_name = query_name.upper()
        self.subject_name = subject_name.upper()

        # Dali info -- update_dali_info
        self.Z_score = None 
        self.rmsd = None
        self.lali = None
        self.nres = None
        self.pID = None
        self.description = None

        # Alignment -- update_aligned_residues
        self.aligned_residues = []

        # Catalytic triad -- update_catalytic_residues
        self.catalytic_nuc = None
        self.catalytic_his = None
        self.catalytic_acid = None

        # Loops -- update_loop_maps
        self.loop_maps = {i: None for i in structure_map if isinstance(i, int)}

        if auto:
            self.auto_calculate(dali_file, align_file, query_pose, subject_pose, 
                catalytic_residues, structure_map, report=report, verbose=verbose)

    def auto_calculate(self, dali_file, align_file, query_pose, subject_pose, 
        catalytic_residues, structure_map, report=False, verbose=False):
        """ Runs all update methods """
        if report:
            lr = '{:<50}{:>50}\n'
            report.write('HEADER: ' + self.subject_name + '\n')
            report.write(lr.format('Query Name:', self.query_name))
            report.write(lr.format('Subject Name:', self.subject_name))
            report.write('\n')
            
        self.update_dali_info(dali_file)
        if report:
            report.write(lr.format('Z score:', self.Z_score))
            report.write(lr.format('RMSD:', self.rmsd))
            report.write(lr.format('LALI:', self.lali))
            report.write(lr.format('NRES:', self.nres))
            report.write(lr.format('pID:', self.pID))
            report.write(lr.format('Description:', self.description))
            report.write('\n')

        self.update_aligned_residues(align_file, query_pose, subject_pose, 
                                     report=report, verbose=verbose)

        self.update_catalytic_residues(catalytic_residues)
        if report:
            res_line = '{:<8}' * 12 + '\n'
            report.write('\n')
            report.write('Catalytic residues\n')
            report.write(res_line.format('Equal', 'Aligned', 
                                         'AA_type', 'PDB_no', 'Pose_no', 'DSSP', 'Bfactor', 
                                         'AA_type', 'PDB_no', 'Pose_no', 'DSSP', 'Bfactor'))
            report.write('Nucleophile: \n')
            if self.catalytic_nuc:
                report.write(self.catalytic_nuc.get_attributes_printout())
            else:
                report.write('None')
            report.write('Histidine: \n')
            if self.catalytic_his:
                report.write(self.catalytic_his.get_attributes_printout())
            else:
                report.write('None')
            report.write('Acid: \n')
            if self.catalytic_acid:
                report.write(self.catalytic_acid.get_attributes_printout())
            else:
                report.write('None')
            report.write('\n')
            report.write('Loops\n')

        self.update_loop_maps(structure_map, query_pose, subject_pose, report=report, verbose=verbose)
        if report:
            report.write('\n')
            report.write('END SUBJECT ' + self.subject_name + '\n')
            report.write('\n'*2)

        return

    def update_dali_info(self, dali_file):
        """
        Update Z_score, rmsd, lali, nres, pID, and description attributes, based 
        on a given dali data file.
        """
        dali_info = get_dali_info(dali_file, self.subject_name)
        self.Z_score = dali_info['Z_score'] 
        self.rmsd = dali_info['rmsd']
        self.lali = dali_info['lali']
        self.nres = dali_info['nres']
        self.pID = dali_info['pID']
        self.description = dali_info['description']
        return

    def update_aligned_residues(self, align_file, query_pose, subject_pose, 
                                report=False, verbose=False):
        """
        Updates aligned_residues property with a list populated from a Dali 
        alignmant file and query and subject poses. The latter are necessary 
        to determine pose numbers and B-factors.
        """
        self.aligned_residues = map_aligned_residues(align_file, 
            self.subject_name, query_pose, subject_pose, 
            report=report, verbose=verbose)
        return

    def update_catalytic_residues(self, catalytic_residues, mode='pdb'):
        """
        Updates the six properties related to identifying the subject protein's 
        catalytic residues, based on a dict of the catalytic residues in the 
        query, in the form of {H: histidine, A: acid residue, N: nucleophile}.
        Default behavior expects that the input dict identifies the catalytic 
        residues by PDB number, not pose number, though this can be altered.
        """
        # Make sure mode is either pdb or pose
        assert mode.lower() in ['pdb', 'pose']

        # Identify he catalytic residues in the subject
        catalytic_residues = map_cat_res(self.aligned_residues, 
            catalytic_residues, mode=mode)

        # Set attribute values
        if catalytic_residues['N']:
            if catalytic_residues['N'].subject_res_type in ['A', 'C', 'S']:
                self.catalytic_nuc = catalytic_residues['N']

        if catalytic_residues['H']:
            if catalytic_residues['H'].subject_res_type in ['H']:
                self.catalytic_his = catalytic_residues['H']
        
        if catalytic_residues['A']:
            if catalytic_residues['A'].subject_res_type in ['D', 'E']:
                self.catalytic_acid = catalytic_residues['A']

        return

    def update_loop_maps(self, structure_map, query_pose, subject_pose, report=False, verbose=False):
        """
        Updates loop_maps attribute with a list of matched_loop objects, based 
        on a structure map and query and subject poses. The structure map should 
        be a dict of the form {'N': range(<N-term_start>, <N-term_end>), 
        1:range(<loop_1_satrt>,<loop_1_end>), ..., 'C': range(<C-term_start>, 
        <C-term_end>)}. The ranges are assumed to be in PDB numbers.
        """
        self.loop_maps = map_structure_elements(structure_map, query_pose, 
            subject_pose, self.aligned_residues, self.subject_name, report=report, verbose=verbose)
        return


def get_dali_info(dali_file, subject):
    """
    Read in appropriate summary from Dali download about this protein, including 
    Z score (indicating structural similarity to the query structure), RMSD to 
    TEV protease (the original query), lali (the number of structurally 
    equivalent CA atoms), nres (the total number of residues in the chain), pID 
    (percentage of identical amino acids in equivalent residues), and PDB 
    description.

    Header line:
    Chain   Z   rmsd lali nres  %id Description
    """
    # Add the hyphen for PDBid-Chain if subject is in form PDBidChain
    if len(subject) == 5:
        subject = '-'.join([subject[:4], subject[-1]])
    
    # Read in Dali summary
    with open(dali_file, 'r') as r:
        match_summaries = r.readlines()

    # Find appropriate line in the summary by PDB name, stopping when found
    summary_line = None
    for ms in match_summaries:
        if subject in ms.upper():
            summary_line = ms.split()
            break

    # Initialize output dict
    dali_info = {}
    dali_info['Z_score'] = None
    dali_info['rmsd'] = None
    dali_info['lali'] = None
    dali_info['nres'] = None
    dali_info['pID'] = None
    dali_info['description'] = None
    
    # If no appropriate line is found, print error message and exit
    if summary_line == None:
        print("No matching protein identified in Dali summary")

    # If line was found, read in its values
    else: 
        dali_info['Z_score'] = summary_line[2]
        dali_info['rmsd'] = summary_line[3]
        dali_info['lali'] = summary_line[4]
        dali_info['nres'] = summary_line[5]
        dali_info['pID'] = summary_line[6]
        dali_info['description'] = ' '.join(summary_line[7:])

    return dali_info


def extract_alignment_block(align_file, subject):
    """
    Reads in a text file with many Dali alignments, and extracts the relevant 
    lines for a given subject. Each alignment block starts with a line that 
    includes 'Z-score', and that string does not appear anywhere else in the 
    blocks but the start line, hence that is used to delineate blocks. 
    """
    # Read in sequence alignment
    with open(align_file, 'r') as r:
        seq_aligns = r.readlines()

    # Find appropriate lines in the summary by PDB name
    begin_block = None
    end_block = None

    for n, sa in enumerate(seq_aligns):
        # Only check starting lines
        if 'Z-score' in sa:
            # Stop capture at the next block after finding the start
            if begin_block != None:
                end_block = n - 1 
                break
    
            # Find beginning of block, where start line includes name
            if subject.upper() in sa.upper():
                begin_block = n
    
    # Extracting relevant text block
    alignment_block = seq_aligns[begin_block:end_block]

    return alignment_block


def get_dali_alignment(align_file, subject):
    """
    Read in sequence alignment file as a set of contiguous strings
    Hacky--may need tweaking to generalize.

    Alignment file has sets of five lines, with each set covering 60 residues 
    in the alignment. 
    The first line (0) is the secondary structure of the query (TEV). 
    The second (1) is the query sequence. 
    The third (2) is the identity match (indicated as positive by |). 
    The fourth (3) is the subject sequence. 
    The fifth (4) is the subject secondary structure.

    Returns these as a dict
    """
    # Extracting approrpiate alignment block from alignment file
    alignment_block = extract_alignment_block(align_file, subject)

    # Cleaning block: delete the first two lines of the block, which are not 
    # alignment, and all blank lines
    abclean = [i.strip() for i in alignment_block[2:] if i != '\n']

    # Chack that there are the right number of lines (a multiple of 5)
    assert len(abclean) % 5 == 0 

    # Concatenating data portions of each alignment line. See docstring.
    align_lines = {0: '', 1: '', 2: '', 3: '', 4: ''}
    for n, line in enumerate(abclean):
        which_set = n % 5
        # Cut off before residue numbers
        if which_set == 0:
            max_len = len(line)
        # Pad short lines
        line_info = line[6:max_len]
        while len(line_info) < max_len - 6:
            line_info += ' '
        # Adding to appropriate set
        align_lines[which_set] += line_info

    # Verifying all lines are equal length, set dict value
    line_lengths = [len(i) for i in align_lines.values()]
    assert all([elem == line_lengths[0] for elem in line_lengths])
    align_lines['length'] = line_lengths[0]

    # Correct dict keys
    align_lines['query_secstruct'] = align_lines.pop(0)
    align_lines['query_sequence'] = align_lines.pop(1)
    align_lines['identity'] = align_lines.pop(2)
    align_lines['subject_sequence'] = align_lines.pop(3)
    align_lines['subject_secstruct'] = align_lines.pop(4)

    return align_lines


def get_posenums_from_dali_str(pose, dali_string, report=False):
    """
    For a given pose and matching string from a Dali alignment, returns a list 
    of pose numbers. Where there ar gaps, the list will include a None entry. 
    Uses Biopython's global alignment with maximum similarity function, with 
    scores from https://towardsdatascience.com 
    Matching characters: +2
    Mismatching character: -1
    Opening a gap: -0.5
    Extending a gap: -0.1 
    """ 
    # If pose has multiple chains, only take first one
    pose_chain_1 = pose.split_by_chain()[1]

    # Get sequence strings from pose and Dali
    ps = pose_chain_1.sequence().upper()    # Pose sequence
    ds = dali_string.upper()        # Dali sequence

    if report:
        report.write('Aligning for pose numbering:\n')
        report.write('Original string:\n')
        report.write(ds)
        report.write('\nPose sequence:\n')
        report.write(ps)
        report.write('\n\n')
        
    # Aligning pose sequence and Dali sequence
    alignments = pairwise2.align.globalms(ps, ds, 2, -1, -0.5, -0.1)

    # Verify that there is only one best alignment
    assert len(alignments) == 1

    # Initializing pose numbering and empty list
    posnum = 1
    pnlist = []
    
    # Filling in list; the first element in the alignment is the pose sequence
    # with dashes inserted to align it with the Dali string
    for a in alignments[0][0]:
        if a == '-':
            pnlist.append(None)
        else:
            assert pose_chain_1.residue(posnum).name1() == a # Verify correct residue 
            pnlist.append(posnum)
            posnum += 1

    # Warning for weird cases like 1CU1, where residues int he pose are placed 
    # differently from the alignment in Dali, possibly due to a reordering or 
    # circular permutation, or a detached chain
    if len(pnlist) > len(dali_string):
        print("Warning: Pose includes residues beyond Dali alignment.")

    # For other edge case, if pose had more residues on the N-term side than 
    # the Dali alignment. If this were to happen, pose numbers would be 
    # incorrect. Break if that happens. The result of such an input would be 
    # extra dashes inserted at the beginning of the second element in the 
    # alignment
    assert alignments[0][1][:len(ds)] == ds 

    if report:
        report.write('Alignment:\n')
        for a in alignments[0]:
            report.write(str(a))
            report.write('\n')
        report.write('\n')

    return pnlist[:len(ds)] 


def map_aligned_residues(align_file, subject, query_pose, subject_pose, 
                         report=False, verbose=False):
    """
    Feed alignment data into a list of aligned_residue objects, each with 
    corresponding information about the position in both query and the subject 
    protease being analyzed.
    """
    # Get alignment 
    alignment = get_dali_alignment(align_file, subject)

    if verbose:
        for k,v in alignment.items():
            report.write(str(k))
            report.write('\n')
            report.write(str(v))
            report.write('\n')
    
    # Get pose number lists, add to alignment
    if verbose:
        alignment['query_pose_numbers'] = get_posenums_from_dali_str(query_pose, 
            alignment['query_sequence'], report=report)
        alignment['subject_pose_numbers'] = get_posenums_from_dali_str(subject_pose, 
            alignment['subject_sequence'], report=report)
    else:
        alignment['query_pose_numbers'] = get_posenums_from_dali_str(query_pose, 
            alignment['query_sequence'])
        alignment['subject_pose_numbers'] = get_posenums_from_dali_str(subject_pose, 
            alignment['subject_sequence'])

    # Initialize aligned residues list
    aligned_residues = []
    
    if report:
        res_line = '{:<8}' * 12 + '\n'
        report.write('Residue alignments:\n')
        report.write(' ' * 16 + '{:<40}{:<40}\n'.format('Query', 'Subject'))
        report.write(res_line.format('Equal', 'Aligned', 
                                     'AA_type', 'PDB_no', 'Pose_no', 'DSSP', 'Bfactor', 
                                     'AA_type', 'PDB_no', 'Pose_no', 'DSSP', 'Bfactor'))

    # Loop through each residue in the alignment, adding aligned_residue
    # objects to the list for each
    for i in range(alignment['length']):
        # Check that at least one pose-numbered residue is present
        q_res = alignment['query_pose_numbers'][i]
        s_res =  alignment['subject_pose_numbers'][i]
        if not any([q_res, s_res]):
            continue

        # Populating aligned_residue object
        a_residue = aligned_residue(i, alignment, query_pose, subject_pose, report=report)

        # Adding aligned_residue object to self.aligned_residues
        aligned_residues.append(a_residue)

    return aligned_residues


def map_cat_res(aligned_residues, catalytic_map, mode='pdb'):
    """
    Using the list of aligned residues, identify the residues in the subject 
    pose that match the query. Requires an input for the catalytic triad in the 
    form of a dict, {H: histidine, A: acid residue, N: nucleophile}, using PDB 
    (not pose) numbering by default. Can use pose numbering instead.
    """
    # Make sure mode is either pdb or pose
    assert mode.lower() in ['pdb', 'pose']

    # Check which mode
    mode_pdb = mode.lower() == 'pdb'
    mode_pose = mode.lower() == 'pose' 

    # Initialize list of matched catalytic residues as all None
    catalytic_residues = {'H': None, 'A': None, 'N': None} 

    # Collect list of just the aligned_residue objects for the catalytic 
    # residues, based on the appropriate pdb/pose numbers
    target_res = catalytic_map.values()
    for ar in aligned_residues:
        # Get appropriate residue number
        if mode_pdb:
            res_number = ar.query_pdb_number
        if mode_pose:
            res_number = ar.query_pose_number

        # Check if that residue number is in the catalytic list
        if res_number in target_res:
            # Check that residues are aligned
            if ar.residues_align:
                # Add to apropriate part of dict
                for k, v in catalytic_map.items():
                    if res_number == v:
                        catalytic_residues[k] = ar

    return catalytic_residues


def map_structure_elements(structure_map, query_pose, subject_pose, 
    aligned_residues, subject_name, res_type='pdb', target='query', report=False, verbose=False):
    """
    Assemble a dict of matched_loop objects, based on the given structure map.
    The structure map should be a dict with the keys listing identified loops, 
    including the unstructured terminal regions. The matched_loop generation 
    will look in the regions between loops for potential splicing sites, so 
    this function takes all loop cases and identifies the residue range within  
    which the matched_loop might look for those sites. Default behavior assumes 
    that the structure_map is in PDB numbers (mode) and that the map is of the 
    query (not subject) protein.
    """
    # Make sure mode is either pdb or pose
    assert res_type.lower() in ['pdb', 'pose']

    # Make sure target is either query or subject
    assert target.lower() in ['query', 'subject']

    # Determine number of loops to match, based on map
    first_loop = min([x for x in structure_map.keys() if isinstance(x, int)])
    last_loop = max([x for x in structure_map.keys() if isinstance(x, int)])

    # Initialize collection dict
    loop_maps = {}

    # Generate a matched_loop for each loop identified in the structure_map
    for loop in structure_map:
        # Get boundaries:
        # One residue past the last residue of upstream loop through
        # One residue before the first residue of downstream loop
        if loop == 'N': # Edge case for N-terminal region (not a loop)
            continue
            n_bound = None
            c_bound = structure_map[1][0] - 1 
            # ^ = Everything up to first res of first loop
        elif loop == 'C': # Edge case for C-terminal region (not a loop)
            continue
            n_bound = structure_map[last_loop][-1] + 1 
            # ^ = Everything after last res of last loop
            c_bound = None
        elif loop == first_loop: # Edge case for N-terminal loop
            n_bound = structure_map['N'][-1] + 1
            c_bound = structure_map[loop + 1][0] - 1
        elif loop == last_loop: # Edge case for C-terminal loop
            n_bound = structure_map[loop - 1][-1] + 1
            c_bound = structure_map['C'][0] - 1
        else: # General case for all interrior loops
            n_bound = structure_map[loop - 1][-1] + 1
            c_bound = structure_map[loop + 1][0] - 1

        # Get subset of aligned_residues between determined bounds
        ar_subset = partition_aligned_residues_list(aligned_residues, 
            n_bound, c_bound, res_type=res_type, target=target)[1]

        # Make matched_loop object and add it to the dict
        #try:
        #    loop_map = matched_loop(subject_name, loop, ar_subset, 
        #        structure_map[loop], query_pose, subject_pose, report=report, verbose=verbose)
        #    loop_maps[loop] = loop_map
        #except Exception as e:
        #    print('{} loop {} failed\n'.format(subject_name, loop))
        #    print(e)
        #    loop_maps[loop] = None
        loop_map = matched_loop(subject_name, loop, ar_subset, 
            structure_map[loop], query_pose, subject_pose, report=report, verbose=verbose)
        loop_maps[loop] = loop_map

    return loop_maps


# Functions used by aligned_residue class or on lists of them

class aligned_residue():
    """
    Data storage structure for a single residue. Includes information about 
    both the target residue in the subject protein and the corresponding aligned 
    residue in the query structure. Information includes secondary structure, 
    whether residues are structurally matched (as opposed to unaligned), and 
    whether the residues are identical. Also stores residue numbers (both PDB
    and pose) for both residues, and B-factor info if available.
    """
    def __init__(self, index, alignment, query_pose, subject_pose, auto=True, report=False):
        # Query attributes 
        self.query_res_type = None
        self.query_sec_struct = None
        self.query_pdb_number = None
        self.query_pose_number = None
        self.query_b_factor = None

        # Subject attributes
        self.subject_res_type = None
        self.subject_sec_struct = None
        self.subject_pdb_number = None
        self.subject_pose_number = None
        self.subject_b_factor = None

        # Determine whether residues are structurally aligned
        self.residues_align = False

        # Determine whether residues are identical
        self.residues_equal = False

        if auto:
            self.auto_populate(index, alignment, query_pose, subject_pose, report=report)

    def auto_populate(self, index, alignment, query_pose, subject_pose, report=False):
        """
        Populates all attributes of the aligned_residue object, based on a 
        given index in an alignment (a dict including strings from a Dali 
        alignment including secondary structure and sequence for both query and 
        subject proteins, and whether the residues are identical; the dict also 
        includes lists of PDB numbers for query and subject--see 
        get_dali_alignment and map_aligned_residues), and query and subject 
        poses.
        """
        # Generate dicts of residue properties for subject and query
        query_res = get_numbers_dss_bfac(index, 
            alignment['query_pose_numbers'], alignment['query_sequence'], 
            alignment['query_secstruct'], query_pose)
        subject_res = get_numbers_dss_bfac(index, 
            alignment['subject_pose_numbers'], alignment['subject_sequence'], 
            alignment['subject_secstruct'], subject_pose)

        # Collect residue identity
        res_identity = alignment['identity'][index]

        # Determine if both query and subjectresidues are present
        query_present = bool(query_res['pose_num'])
        subject_present = bool(subject_res['pose_num'])

        # Determine whether residues are structurally aligned, based on case
        if all([query_present, subject_present]):
            if   all([i['sequence'].isupper() for i in [query_res, subject_res]]):
                residues_align = True
            elif all([i['sequence'].islower() for i in [query_res, subject_res]]):
                residues_align = False
            else:
                print('Residue cases do not match')
                print(query_res['pose_num'], query_res['sequence'], 
                    subject_res['pose_num'], subject_res['sequence'])
                assert False
        else:
            residues_align = False

        # Determine res identity, based on whether connection line was drawn
        if res_identity == '|':
            residues_equal = True
            assert query_res['sequence'] == subject_res['sequence']
        else:
            residues_equal = False

        # Update attributes
        if query_present:
            self.query_res_type = query_res['sequence'].upper()
            self.query_sec_struct = query_res['secstruct']
            self.query_pdb_number = query_res['pdb_num']
            self.query_pose_number = query_res['pose_num']
            self.query_b_factor = query_res['b_factor']

        if subject_present:
            self.subject_res_type = subject_res['sequence'].upper()
            self.subject_sec_struct = subject_res['secstruct']
            self.subject_pdb_number = subject_res['pdb_num']
            self.subject_pose_number = subject_res['pose_num']
            self.subject_b_factor = subject_res['b_factor']

        self.residues_align = residues_align
        self.residues_equal = residues_equal
        
        # Report
        if report:
            report.write(self.get_attributes_printout())

        return

    def get_attributes_printout(self):
        """ 
        Makes a single string will all residue attributes. 
        Column 1:  Whether the residues are equivalent
        Column 2:  Whether the residues are aligned
        Column 3:  Query amino acid type
        Column 4:  Query PDB number
        Column 5:  Query Pose number
        Column 6:  Query secondary structure
        Column 7:  Query B factor
        Column 8:  Subject amino acid type
        Column 9:  Subject PDB number
        Column 10: Subject Pose number
        Column 11: Subject secondary structure
        Column 12: Subject B factor
        """
        # Template
        res_line = '{:<8}' * 12 + '\n'

        # Round B factors, if they are numerical
        if self.query_b_factor:
            qb = round(self.query_b_factor, 3)
        else:
            qb = None
        if self.subject_b_factor:
            sb = round(self.subject_b_factor, 3)
        else:
            sb = None

        # List attributes
        attribs = [bool(self.residues_equal), bool(self.residues_align),
            self.query_res_type, self.query_pdb_number, 
            self.query_pose_number, self.query_sec_struct, qb,
            self.subject_res_type, self.subject_pdb_number, 
            self.subject_pose_number,self.subject_sec_struct, sb]

        # Convert all attributes to a string
        out_as_str = [str(i) for i in attribs]

        # Return formatted output
        return res_line.format(*out_as_str)


def get_numbers_dss_bfac(index, pose_numbers, seq_string, secstruct_string, 
    pose):
    """
    For a given index, list of pose numbers, protein sequence string, secondary  
    structure string, and pose, pull the secondary structure and sequence at the 
    given index and determine the PDB and pose numbers and B factor from the pose. It is 
    possible that there are gaps in the strings, which correspond to nothing in 
    the pose. In these cases, will return None for all values. Returns a dict 
    with these five values. 
    """
    # Initialize dict
    res_specs = {'sequence': None, 'secstruct': None, 'pdb_num': None, 
        'pose_num': None, 'b_factor': None}

    # Check whether there is a pose number for this residue
    pose_number = pose_numbers[index]

    # If there is, update values
    if pose_number:
        # Residue letter--case left unchanged to check alignment
        res_specs['sequence'] = seq_string[index]

        # Secondary structure
        res_specs['secstruct'] = secstruct_string[index].upper()
        
        # PDB number
        pdb_number = pose.pdb_info().pose2pdb(pose_number).split()[0]
        res_specs['pdb_num'] = int(pdb_number)

        # Pose number
        res_specs['pose_num'] = pose_number

        # Get query B-factor
        res_specs['b_factor'] = get_b_factor(pose, pose_number)

    return res_specs


def get_ar_property(aligned_residue, target, attribute):
    """
    Converts casual user input into a request for a specific property of an 
    aligned_res object for the get_res function. Target is either query (q) or
    subject (s), and is case insensitive. Attribute is residue type (rt, type, 
    res_type), secondary structure (ss, secstruct, sec_struct), pdb number (pd, 
    pdb, pdbnum, pdb_num, pdb_number), pose number (po, pos, pose, posenum, 
    pose_num, pose_number), or B factor (b, bfac, b_factor), and is case 
    insensitive.
    """
    # List acceptable inputs for target
    targets = {'query': ['Q', 'QUER', 'QUERY'], 
               'subject': ['S', 'SUBJ', 'SUBJECT']}

    # List acceptable inputs for attribute 
    attributes = {'res_type':    ['RT', 'TYPE', 'RES_TYPE', 'RESIDUE_TYPE'],
                  'sec_struct':  ['SS', 'SECSTRUCT', 'SEC_STRUCT', 
                                  'SECONDARY_STRUCTURE'],
                  'pdb_number':  ['PD', 'PDB', 'PDBNUM', 'PDB_NUM', 
                                  'PDB_NUMBER'],
                  'pose_number': ['PO', 'POS', 'POSE', 'POSENUM', 
                                  'POSE_NUMBER'],
                  'b_factor':    ['B', 'BFAC', 'B_FACTOR']}

    # Determine target residue type
    target = target.upper()
    out_target = None
    for t, v in targets.items():
        if target in v:
            out_target = t
            break

    # Determine attribute type correct,
    # replacing dashes and spacees with underscores to accomodate inputs
    attribute = attribute.upper().replace('-','_').replace(' ', '_')
    out_attribute = None
    for a, v in attributes.items():
        if attribute in v:
            out_attribute = a
            break

    # Confirm that the inputs produced a valid target and attribute
    invalid = False
    if not out_target:
        print("Invalid input for target") 
        invalid = True
    if not out_attribute:
        print("Invalid input for attribute") 
        invalid = True
    if invalid:
        assert False

    retrieve = '_'.join([out_target, out_attribute])
    return getattr(aligned_residue, retrieve)


def find_ar(seek_value, aligned_residues, target, attribute):
    """

    """
    # List acceptable inputs for target
    targets = {'query': ['Q', 'QUER', 'QUERY'], 
               'subject': ['S', 'SUBJ', 'SUBJECT']}

    # List acceptable inputs for attribute 
    attributes = {'res_type':    ['RT', 'TYPE', 'RES_TYPE', 'RESIDUE_TYPE'],
                  'sec_struct':  ['SS', 'SECSTRUCT', 'SEC_STRUCT', 
                                  'SECONDARY_STRUCTURE'],
                  'pdb_number':  ['PD', 'PDB', 'PDBNUM', 'PDB_NUM', 
                                  'PDB_NUMBER'],
                  'pose_number': ['PO', 'POS', 'POSE', 'POSENUM', 
                                  'POSE_NUMBER'],
                  'b_factor':    ['B', 'BFAC', 'B_FACTOR']}

    # Determine target residue type
    target = target.upper()
    out_target = None
    for t, v in targets.items():
        if target in v:
            out_target = t
            break

    # Determine attribute type correct,
    # replacing dashes and spacees with underscores to accomodate inputs
    attribute = attribute.upper().replace('-','_').replace(' ', '_')
    out_attribute = None
    for a, v in attributes.items():
        if attribute in v:
            out_attribute = a
            break

    # Confirm that the inputs produced a valid target and attribute
    invalid = False
    if not out_target:
        print("Invalid input for target") 
        invalid = True
    if not out_attribute:
        print("Invalid input for attribute")
    if invalid:
        assert False

    retrieve = '_'.join([out_target, out_attribute])
    for ar in aligned_residues:
        if getattr(ar, retrieve) == seek_value:
            return ar


def partition_aligned_residues_list(aligned_residues, n_bound, c_bound, 
    target='query', res_type='pdb'):
    """
    Returns a subsets of a list of aligned_residue objects, partitioned with a  
    given set of boundaries. The boundaries can be either PDB or pose numbers.  
    The bounds can correspond to residues in either the query protein or the 
    subject. Returns a list of three lists: the set before, the set between, 
    and the set after. The bounds are included in the middle set.
    """
    # Loop through list to find boundary residues for given res_type and target
    for n, ar in enumerate(aligned_residues):
        resnum = get_ar_property(ar, target, res_type)
        if resnum == n_bound:
            start_index = n
        if resnum == c_bound:
            end_index = n

    # Address edge cases where no boundary is given, overwriting 
    if n_bound == None:
        start_index = 0

    if c_bound == None: 
        end_index = len(aligned_residues)

    # Define lists
    set_before = aligned_residues[:start_index]
    middle_set = aligned_residues[start_index : end_index + 1]
    set_after = aligned_residues[end_index + 1:]

    return [set_before, middle_set, set_after]


def find_nearest_matches(aligned_residues):
    """
    Scans through a list of aligned_residue objects and identifies, starting 
    from the beginning of the list, the first matching residue, the first 
    residue that matches and is also a B-sheet, the last residue that matches, 
    and the last residue that matches and is a B-sheet. Returns these as a dict.
    It is possible that no matches will be found, in which case the dict will 
    include None.
    """
    # Identify aligned residues
    align_mask = [ar.residues_align for ar in aligned_residues]
    aligning_elements = [n for n, i in enumerate(align_mask) if i]

    if len(aligning_elements) == 0:
        return None, None

    # Trim list to avoid any unaligned residues
    first_align = aligning_elements[0]
    last_align = aligning_elements[0]
    if len(aligning_elements) > 1:
        for n, i in enumerate(aligning_elements):
            if n == 0:
                continue
            if i == aligning_elements[n - 1] + 1:
                last_align = i
            else:
                break

    return aligned_residues[first_align], aligned_residues[last_align] 


def get_residue_list(aligned_residues, target='query', res_type='pdb'):
    """
    For a given list of aligned_residue objects, returns a list of residues of 
    the specified type. Default behavior is to return the PDB numbers of the 
    query protein, though the subject protein can be searched instead, and the 
    pose numbers or B factors can be searched instead.
    """
    # Initialize output list
    out_res_list = []

    # Loop through list and collect the desired attribute
    for ar in aligned_residues:
        res_retrieve = get_ar_property(ar, target, res_type)
        out_res_list.append(res_retrieve)

    return out_res_list


def prealign_by_dali(query_pose, subject_pose, aligned_residues):
    """
    For a subject pose and a query pose, reads the list of aligned_residues and 
    creates selectors for both poses of only the residues that Dali identified 
    as aligned, then aligns based on those residues.
    """
    # Make list of only the residues that are aligned
    align_mask = [res.residues_align for res in aligned_residues]
    only_aligned_residues = apply_mask(aligned_residues, align_mask)
            
    # Get pose numbers for query and subject
    query_residues = get_residue_list(only_aligned_residues, 
                                      res_type='pose', target='query')
    subject_residues = get_residue_list(only_aligned_residues, 
                                        res_type='pose', target='subject')
    
    # Align poses
    rmsd = align_protein_sections(query_pose, query_residues, 
                           subject_pose, subject_residues, mode='BB')
    
    return rmsd


# Functions used by matched_loop class

class matched_loop():
    """
    Data storage structure for loops. When taking in a loop of a query 
    structure, finds the edges bordering it (usually B-sheets) and looks for 
    residue matches within the given boundaries, which should be the starts of
    the next loops. Input residues should use PDB (not pose) numbers.
    """
    def __init__(self, source, l_name, local_residues, l_range, query_pose, 
        subject_pose, auto=True, report=False, verbose=False, dump=False):

        # Loop source and name
        self.loop_source = source
        self.loop_name = l_name

        # Flanking residues -- update_matches
        self.n_farthest_match = None
        self.n_nearest_match = None
        self.c_nearest_match = None
        self.c_farthest_match = None
        self.flanking_matches_found = False
        self.query_loop_size = None
        self.subject_loop_size = None
        self.loop_size_change = None

        # Best matched residues for loop swap -- update_splice_sites
        #self.N_outside_overlap_residue = None
        #self.N_splice_residue = None
        #self.N_loop_end_residue = None
        #self.C_loop_end_residue = None
        #self.C_splice_residue = None
        #self.C_outside_overlap_residue = None
            # Numerical properties of loop
        #self.N_overlap_size = None
        #self.C_overlap_size = None

        # Loop RMSD is query and subject are the same size -- update_loop_rmsd
        self.rmsd = None

        # Average B factor for the loop -- update_b_factor
        self.b_factor = None

        # Check for discontinuities -- update_continuity
        self.is_continuous = None

        # Loop proximity to peptide substrate -- update_proximity_to_substrate
        #self.closest_residue_distance = None
        #self.close_substrate_residues = []
        
        # Loop alignment ensemble -- update_alignment_ensemble
        self.alignment_ensemble = []

        # Evaluate whether loop is a suitable target -- update_suitability
        self.is_n_match = None
        self.is_c_match = None
        self.is_near_target = None
        self.is_not_domain = None
        self.is_different_from_original = None
        self.is_ordered = None
        self.is_possible_target = None
        
        if auto:
            self.auto_calculate(local_residues, l_range, query_pose, subject_pose, 
                                report=report, verbose=verbose, dump=dump)

    def auto_calculate(self, local_residues, l_range, query_pose, subject_pose, 
        report=False, verbose=False, dump=False):
        """ Runs all update methods """
        if report:
            lr = '{:<50}{:>50}\n'
            res_line = '{:<8}' * 12 + '\n'
            report.write(lr.format('Source PDB:', self.loop_source))
            report.write(lr.format('Loop site:', self.loop_name))
            report.write('\n')
            report.write('Local residues:\n')
            report.write(' ' * 16 + '{:<40}{:<40}\n'.format('Query', 'Subject'))
            report.write(res_line.format('Equal', 'Aligned', 
                                         'AA_type', 'PDB_no', 'Pose_no', 'DSSP', 'Bfactor', 
                                         'AA_type', 'PDB_no', 'Pose_no', 'DSSP', 'Bfactor'))
            for res in local_residues:
                report.write(res.get_attributes_printout())
            report.write('\n')
            
        res_sets = self.update_matches(local_residues, l_range[0], l_range[-1])
        if report:
            report.write(lr.format('Flanking matches found:', str(self.flanking_matches_found)))
            report.write('Farthest N-terminal match:\n')
            if self.n_farthest_match:
                report.write(self.n_farthest_match.get_attributes_printout())
            else:
                report.write('None\n')
            report.write('Nearest N-terminal match:\n')
            if self.n_nearest_match:
                report.write(self.n_nearest_match.get_attributes_printout())
            else:
                report.write('None\n')
            report.write('Nearest C-terminal match:\n')
            if self.c_nearest_match:
                report.write(self.c_nearest_match.get_attributes_printout())
            else:
                report.write('None\n')
            report.write('Farthest C-terminal match:\n')
            if self.c_farthest_match:
                report.write(self.c_farthest_match.get_attributes_printout())
            else:
                report.write('None\n')
            report.write('\n')

        if report:
            if not self.flanking_matches_found:
                report.write('Skipped alignment window check due to lack of matched residues.\n')
                report.write('Skipped RMSD comparison due to lack of matched residues.\n')
                report.write('Skipped B factor calculation due to lack of matched residues.\n')
            else:    
                report.write('Alignment window check:\n')
        #n_splice_range, c_splice_range = self.update_splice_sites(res_sets, query_pose, subject_pose, 
        #                                                          report=report, verbose=verbose, dump=dump)
        if report:
            if self.flanking_matches_found:
                if self.loop_size_change != 0:
                    report.write('Skipped RMSD comparison since subject loop is not the same size as query.\n')
                    report.write('\n')
                
        self.update_loop_rmsd(query_pose, subject_pose)
        self.update_b_factor(local_residues)
        self.update_continuity(subject_pose, report=report)

        #if all([n_splice_range, c_splice_range]):
        #    self.update_proximity_to_substrate(query_pose, subject_pose, 
        #                                       n_splice_range + c_splice_range)

        if report:
            ht = '{:<6}' * 9 + '{:<12}' + '{:<6}' * 2 + '\n'
            lc_head = ['Nout', 'Nsplc', 'Csplc', 'Cout', 'novrlp', 'covrlp', 
                      'rmsd', 'clashs', 'contct', 'subst_res', 'near', 'fit']
            report.write(ht.format(*lc_head))
        self.update_alignment_ensemble(res_sets, query_pose, subject_pose, report=report)
        
        self.update_suitability()
        if report:
            report.write(lr.format('\nLOOP IS CANDIDATE FOR EXCHANGE:', self.is_possible_target))
            report.write('\n')
            report.write('N splice site:\n')
            if self.n_nearest_match:
                report.write( self.n_nearest_match.get_attributes_printout())
            else:
                report.write('None\n')
            report.write('C splice site:\n')
            if self.c_nearest_match:
                report.write( self.c_nearest_match.get_attributes_printout())
            else:
                report.write('None\n')
            #report.write(lr.format('N-side overlap length:', str(self.N_overlap_size)))
            #report.write(lr.format('C-side overlap length:', str(self.C_overlap_size)))
            report.write('\n')
            #if self.closest_residue_distance:
            #    crd = round(self.closest_residue_distance, 3)
            #else:
            #    crd = 'None'
            #report.write(lr.format('Closest CA distance to substrate:', crd))
            #if self.close_substrate_residues:
            #    cr_list = ','.join([str(i) for i in self.close_substrate_residues])
            #else:
            #    cr_list = 'None'
            #report.write(lr.format('Substrate residues within 8A of loop:', cr_list))
            report.write(lr.format('Loop has residues within 9A of substrate:', self.is_near_target))
            report.write(lr.format('Loop does not clash:', self.is_not_clash))
            report.write('\n')
            report.write(lr.format('Query loop length:', str(self.query_loop_size)))
            report.write(lr.format('Subject loop length:', str(self.subject_loop_size)))
            report.write(lr.format('Loop is not a full domain (>50):', self.is_not_domain))
            report.write('\n')
            report.write(lr.format('Difference in subject vs query length:', str(self.loop_size_change)))
            report.write(lr.format('Loop RMSD:', str(self.rmsd)))
            report.write(lr.format('Loop is structurally asimilar:', self.is_different_from_original))
            report.write('\n')
            if self.b_factor: 
                bf = round(self.b_factor, 3)
            else:
                bf = 'None'
            report.write(lr.format('Loop average B factor:', bf))
            if self.is_ordered: 
                ordered = self.is_ordered
            else:
                ordered = 'None'
            report.write(lr.format('Average B factor is low:', ordered))
            report.write('\n')
            if self.is_continuous: 
                continuous = self.is_continuous
            else:
                continuous = 'None'
            report.write(lr.format('Loop is not missing residues:', continuous))
            report.write('\n'*2)

        return

    def update_matches(self, local_residues, loop_start, loop_end):
        """
        Update n_ and c_ nearest_match and properties, identifying the closest 
        and farthest matching residues flanking the loop. Inputs include a 
        list of aligned_residue objects that should be trimmed down from the 
        full alignment to just the region between the adjacent loops, and the 
        PDB numbers of the start and end of the target loop, which should be 
        within the listed residues.
        """
        # Partition local residue set at loop start and end points
        res_sets = partition_aligned_residues_list(local_residues,  
            loop_start, loop_end, res_type='pdb', target='query')

        # Collect matching residues, reversing the C-side set so as to count 
        # toward the loop.
        n_far, n_near = find_nearest_matches(res_sets[0])
        c_far, c_near = find_nearest_matches(res_sets[2][::-1])

        # Update properties 
        self.n_farthest_match = n_far
        self.n_nearest_match = n_near
        self.c_nearest_match = c_near
        self.c_farthest_match = c_far

        # Indicate if matches were found
        if self.n_nearest_match and self.c_nearest_match:
            self.flanking_matches_found = True
            
            # Determine loop sizes
            self.query_loop_size = get_ar_property(self.c_nearest_match, 'q', 'po') -                 get_ar_property(self.n_nearest_match, 'q', 'po') + 1
            self.subject_loop_size = get_ar_property(self.c_nearest_match, 's', 'po') -                 get_ar_property(self.n_nearest_match, 's', 'po') + 1
            self.loop_size_change = self.subject_loop_size - self.query_loop_size

        return res_sets
        
    def update_splice_sites(self, res_sets, query_pose, subject_pose, 
                            report=False, verbose=False, dump=False, out_name=None):
        """
        Updates splice residues, where residues from the query protein should be 
        replaced with residues from the subject protein. If no overlapping 
        residues have been found for this loop between the query and subject, 
        splice sites will remain None. If matches have been found, the ranges of 
        matched residues (both on the N-side and C-side) will be compared in 
        sliding windows to determine the set of residues that yield the lowest 
        backbone RMSD between the subject and query. The set that yields the 
        lowest RMSD will be taken as the splice residues. Updates the splice 
        N_ and C_ splice_residue and overlap_size attributes, as well as the 
        query_ and subject_loop_size attributes.
        """
        # Nothing to update if no matching residues were identified
        if not self.flanking_matches_found:
            return None, None        

        # Generate sub-lists for the overlap windows
        n_start = res_sets[0].index(self.n_farthest_match)
        n_end = res_sets[0].index(self.n_nearest_match) 
        n_window = res_sets[0][n_start: n_end + 1]

        c_start = res_sets[2].index(self.c_nearest_match)
        c_end = res_sets[2].index(self.c_farthest_match) 
        c_window = res_sets[2][c_start: c_end + 1]

        # Find windows of best-fitting residues
        n_splice_range, c_splice_range = find_splice_sites(n_window, c_window, 
            query_pose, subject_pose, report=report, verbose=verbose, dump=dump,
            out_name=out_name)
        
        # Find spliced loop termini
        all_res = res_sets[0] + res_sets[1] + res_sets[2]
        N_splice = n_splice_range[-1]                      # Closest to loop
        C_splice = c_splice_range[0]                       # Closest to loop
        N_loop_end = None
        C_loop_end = None
        for res in all_res:
            if res.subject_pose_number == N_splice.subject_pose_number + 1:
                N_loop_end = res
            if res.subject_pose_number == C_splice.subject_pose_number - 1:
                C_loop_end = res

        # Update attributes 
        self.N_splice_residue = N_splice                    
        self.C_splice_residue = C_splice                    
        self.N_loop_end_residue = N_loop_end
        self.C_loop_end_residue = C_loop_end
        self.N_outside_overlap_residue = n_splice_range[0]  # Farthest from loop
        self.C_outside_overlap_residue = c_splice_range[-1] # Farthest from loop
        self.N_overlap_size = len(n_splice_range)
        self.C_overlap_size = len(c_splice_range)
        self.query_loop_size = self.C_splice_residue.query_pose_number -             self.N_splice_residue.query_pose_number - 1
        self.subject_loop_size = self.C_splice_residue.subject_pose_number -              self.N_splice_residue.subject_pose_number - 1
        self.loop_size_change = self.subject_loop_size - self.query_loop_size

        return n_splice_range, c_splice_range

    def update_loop_rmsd(self, query_pose, subject_pose):
        """ 
        If the subject and query loops are the same size, take CA RMSD 
        """
        # Nothing to update if no matching residues were identified
        if not self.flanking_matches_found:
            return 
        
        # Nothing to update if the loops are not the same size
        if self.loop_size_change != 0:
            return

        # Getting lists of pose numbers for the query and subject loops between 
        # the splice sites
        query_list = list(range(get_ar_property(self.n_nearest_match, 'q', 'po'), 
                           get_ar_property(self.c_nearest_match, 'q', 'po') + 1))
        subject_list = list(range(get_ar_property(self.n_nearest_match, 's', 'po'), 
                             get_ar_property(self.c_nearest_match, 's', 'po') + 1))
        
        # Calculate RMSD
        self.rmsd = align_protein_sections(query_pose, query_list, 
            Pose(subject_pose), subject_list)

        return

    def update_b_factor(self, local_residues):
        """
        Updates the b_factor attribute with the average backbone B factor of 
        all residues in the loop.
        """
        # Nothing to update if no matching residues were identified
        if not self.flanking_matches_found:
            return 
        
        # Get subset of residue range
        first_loop_res = local_residues.index(self.n_nearest_match) + 1
        last_loop_res = local_residues.index(self.c_nearest_match)
        loop_residues = local_residues[first_loop_res:last_loop_res]

        # Collect list of B factors 
        b_factors = get_residue_list(loop_residues, res_type='b', 
            target='subject')

        # Exclude empty values
        b_factors_clean = [i for i in b_factors if isinstance(i, float)]

        # Take average backbone B factor for the loop
        b_average = np.average(b_factors_clean)

        # Update attribute
        self.b_factor = b_average

        return

    def update_continuity(self, subject_pose, report=False):
        """
        Checks whether the loop includes chain breaks. Updates the is_continuous 
        attribute.
        """
        # Nothing to update if no matching residues were identified
        if not self.flanking_matches_found:
            return

        # Create loop subpose
        subpose_start = self.n_nearest_match.subject_pose_number
        subpose_end = self.c_nearest_match.subject_pose_number
        subpose = Pose(subject_pose, subpose_start, subpose_end)

        # Check for breaks
        continuous, c_n_distances, break_sites = check_pose_continuity(subpose)
        
        if report:
            res_in_subpose = range(subpose_start, subpose_end)
            assert len(c_n_distances) == len(res_in_subpose)
            report.write("Continuity check: C-N distances\n")
            for n, i in enumerate(c_n_distances):
                report.write('{} \t {}\n'.format(res_in_subpose[n], round(i, 3)))
            report.write('\n')

        # Announce identified breaks
        if not continuous:
            print('Chain broken:')
            for b in break_sites:
                print('\t', subpose.pdb_info().pose2pdb(b).split()[0])

        # Update attribute
        self.is_continuous = continuous

        return

    def update_alignment_ensemble(self, res_sets, query_pose, subject_pose,
                                 report=None):
        """

        """
        # Nothing to update if no matching residues were identified
        if not self.flanking_matches_found:
            return 
        
        # Abort this (most time-consuming) step if other conditions are no
        if self.subject_loop_size:
            res_count_check = self.subject_loop_size <= 50
        else: 
            res_count_check = True
        if self.loop_size_change == 0:
            similarity_check = self.rmsd > 0.5
        else:
            similarity_check = True
        if self.b_factor:
            b_factor_check = self.b_factor <= 50
        else:
            b_factor_check = True
        if not all([res_count_check, similarity_check, b_factor_check, 
            self.is_continuous]):
            return

        # Generate sub-lists for the overlap windows
        n_start = res_sets[0].index(self.n_farthest_match)
        n_end = res_sets[0].index(self.n_nearest_match) 
        n_window = res_sets[0][n_start: n_end + 1]

        c_start = res_sets[2].index(self.c_nearest_match)
        c_end = res_sets[2].index(self.c_farthest_match) 
        c_window = res_sets[2][c_start: c_end + 1]

        self.alignment_ensemble = make_alignment_ensemble_for_loop(
            n_window, c_window, query_pose, subject_pose, self, report=report)

        return

    def update_proximity_to_substrate(self, query_pose, subject_pose, overlaps):
        """
        Updates the attributes relating to the subject loop's proximity to the 
        query's substrate peptide. This is for purposes of weeding out loops 
        that are not close enough to the substrate to directly affect 
        specificity. Determines the coordinates of the substrate in the query 
        and, after aligning the substrate appropriately with the query, the 
        coordinates of the loop in the subject. Updates loop_near_substrate, 
        closest_residue_distance, and close_substrate_residues.
        """
        # Nothing to update if no matching residues were identified
        if not self.flanking_matches_found:
            return 
        
        # Getting coordinates
        loop_start = self.N_splice_residue.subject_pose_number + 1
        loop_end = self.C_splice_residue.subject_pose_number - 1
        substrate_coords, loop_coords = get_substrate_and_loop_coords(
            query_pose, subject_pose, overlaps, loop_start, loop_end)

        # Determine closest proximity and list of nearby substrate residues
        closest_distance, nearby_substrate_residues =             find_proximity_to_substrate (substrate_coords, loop_coords)

        # Update attributes
        self.closest_residue_distance = closest_distance
        self.close_substrate_residues = nearby_substrate_residues 

        return

    def update_suitability(self, max_res_count=50, min_rmsd=0.5, 
        max_b_factor=50):
        """
        Runs a set of checks to determine whether the loop may be viable for 
        substitution. Rejects loops on the basis of being too far away from the
        substrate to interact, being too large (avoiding domain insertions), 
        being of the same size as the original loop with a small RMSD, lacking 
        matched residues flanking the loop (cannot graft), or having gaps in the 
        crystal structure.
        """
        # Check that there is a matched residue on the N-terminal side
        if self.loop_name == 'N':
            n_match_check = True
        else:
            n_match_check = bool(self.n_nearest_match)

        # Check that there is a matched residue on the C-terminal side
        if self.loop_name == 'C':
            c_match_check = True
        else:
            c_match_check = bool(self.c_nearest_match)

        # Check that loop residue count is not too large
        if self.subject_loop_size:
            res_count_check = self.subject_loop_size <= max_res_count
        else: 
            res_count_check = True

        # Check that if loop is the same size as the query, that it is 
        # structurally different, based on backbone RMSD
        if self.loop_size_change == 0:
            similarity_check = self.rmsd > min_rmsd
        else:
            similarity_check = True

        # B factor check
        if self.b_factor:
            b_factor_check = self.b_factor <= max_b_factor
        else:
            b_factor_check = True

        # Check that there are residues within range of substrate
        plcs_near_target = [c for c in self.alignment_ensemble if c.is_near_target]
        proximity_check = bool(len(plcs_near_target))
        
        # Check that there are non-clashing positions of the loop
        plcs_no_clase = [c for c in self.alignment_ensemble if c.is_low_clash]
        clash_check = bool(len(plcs_no_clase))
        
        # Updating attributes
        self.is_n_match = n_match_check
        self.is_c_match = c_match_check
        self.is_not_domain = res_count_check
        self.is_different_from_original = similarity_check
        self.is_ordered = b_factor_check
        self.is_near_target = proximity_check
        self.is_not_clash = clash_check
        self.is_possible_target = all([n_match_check, c_match_check, 
            proximity_check, res_count_check, similarity_check, b_factor_check, 
            self.is_continuous, clash_check])

        return


def find_splice_sites(n_window, c_window, query_pose, subject_pose, 
                      report=False, verbose=False, dump=False, out_name=None):
    """
    For a given pair of residue windows on the N-side and C-side of a loop, 
    break each window into subframes of length varied between 1 and the full 
    window size. Then take those subframes from each and align based on the 
    residues in those subframes, and determine the best backbone RMSD score.
    The combination of residues that yields the lowest RMSD will be returned.
    """
    if report:
        lf = '{:<45}{:<45}{:>10}\n'
        full_n = get_residue_list(n_window, mode='pose', target='subject')
        full_c = get_residue_list(c_window, mode='pose', target='subject')
        nwin = ','.join([str(i) for i in full_n])
        cwin = ','.join([str(i) for i in full_c])
        report.write(lf.format(nwin, cwin, ''))
        report.write('\n')
        report.write(lf.format('Left Frame', 'Right Frame', 'RMSD'))
        
    # Prioritize B-sheet splice sites
    #is_n_beta = False
    #n_beta = [i for i in n_window if i.subject_sec_struct == 'E']
    #if n_beta:
    #    n_frame = n_beta
    #    is_n_beta = True
    #else:
    #    n_frame = n_window
    #is_n_beta = True
    #n_frame = n_window    
    #is_c_beta = False
    #c_beta = [i for i in c_window if i.subject_sec_struct == 'E']
    #if c_beta:
    #    c_frame = c_beta
    #    is_c_beta = True
    #else:
    #    c_frame = c_window
    #is_c_beta = True
    #c_frame = c_window    
    #if report:
        #report.write('Beta-only windows\n')
        #if is_n_beta:
        #    beta_n = get_residue_list(n_frame, mode='pose', target='subject')
        #    nbwin = ','.join([str(i) for i in beta_n])
        #else:
        #    nbwin = 'None'
        #if is_c_beta:
        #    beta_c = get_residue_list(c_frame, mode='pose', target='subject')
        #    cbwin = ','.join([str(i) for i in beta_c])
        #else:
        #    cbwin = 'None'
        #report.write(lf.format(nbwin, cbwin, ''))
        #report.write('\n')
        #report.write(lf.format('Left Frame', 'Right Frame', 'RMSD'))
        
    # List subframes within overlap windows
    n_subframes = variable_sliding_window(n_window)
    c_subframes = variable_sliding_window(c_window)

    # Initialize collection objects for best alignment
    best_rmsd = 1000 # Arbitrarily high
    best_n_set = n_subframes[-1] # Defaults to the full alignment region
    best_c_set = c_subframes[-1] # Defaults to the full alignment region

    # Loop through all subframes and align based on those residues, taking 
    # the set with the lowest RMSD
    for ns in n_subframes:
        for cs in c_subframes:
            # Do not allow two single-residue frames--minimum three residues to align
            if len(ns) + len(cs) == 2:
                continue
            else: 
                # Get pose residues list for the residues to take RMSD
                q_list_n = get_residue_list(ns, mode='pose', target='query')
                q_list_c = get_residue_list(cs, mode='pose', target='query')
                s_list_n = get_residue_list(ns, mode='pose', target='subject')
                s_list_c = get_residue_list(cs, mode='pose', target='subject')

                # Align these residues and get RMSD
                aligned_subject = Pose(subject_pose)
                rmsd = align_protein_sections(query_pose, q_list_n + q_list_c, 
                    aligned_subject, s_list_n + s_list_c)

                # Compare RMSD to best, update best if lower
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_n_set = ns
                    best_c_set = cs

                if report and verbose:
                    nwin = ','.join([str(i) for i in s_list_n])
                    cwin = ','.join([str(i) for i in s_list_c])
                    report.write(lf.format(nwin, cwin, round(rmsd,5)))

                if dump:
                    subpose = Pose(aligned_subject, s_list_n[0], s_list_c[-1])
                    overlap_bounds = [s_list_n[0], s_list_n[-1], s_list_c[0], s_list_c[-1]]
                    if out_name:
                        outpdb = out_name + '_N{}-{}_C{}-{}.pdb'.format(*[str(i) for i in overlap_bounds])
                    else:
                        outpdb = 'aligned_loop_N{}-{}_C{}-{}.pdb'.format(*[str(i) for i in overlap_bounds])
                    subpose.dump_pdb(outpdb)                
    
    if report:
        best_n_list = get_residue_list(best_n_set, mode='pose', target='subject')
        best_c_list = get_residue_list(best_c_set, mode='pose', target='subject')
        nbest = ','.join([str(i) for i in best_n_list])
        cbest = ','.join([str(i) for i in best_c_list])
        report.write('Best frames:\n')
        report.write(lf.format(nbest, cbest, round(best_rmsd, 3)))
        report.write('\n')

    return best_n_set, best_c_set


def get_substrate_and_loop_coords(query_pose, subject_pose, overlaps, 
    loop_start, loop_end):
    """
    Given a query pose with two chains, where the second chain is the substrate, 
    and a subject pose, the list of overlapping residues to superimpose the 
    loop of the substrate onto the query, and the pose residue numbers that 
    start and end the loop on the subject, returns lists of the CA coordinates 
    of the substrate and the subject loop.
    """
    # Populating list of CA coordinates of substrate peptide
    substrate_chain = query_pose.split_by_chain()[2]
    substrate_coords = list_pose_coords(substrate_chain)

    # Get pose residues list for the overlap residues
    q_list = get_residue_list(overlaps, mode='pose', target='query')
    s_list = get_residue_list(overlaps, mode='pose', target='subject')

    # Creating copy of subject pose and aligning its loop to the query
    aligned_subject = Pose(subject_pose)
    rmsd = align_protein_sections(query_pose, q_list, aligned_subject, s_list)

    # Making loop-only subpose
    aligned_loop = Pose(aligned_subject, loop_start, loop_end)

    # Populating list of CA coordinates of loop
    loop_coords = list_pose_coords(aligned_loop)

    return substrate_coords, loop_coords


def make_alignment_ensemble_for_loop(n_window, c_window, 
    query_pose, subject_pose, matched_loop, report=None):
    """

    """
    # List subframes within overlap windows
    n_subframes = variable_sliding_window(n_window)
    c_subframes = variable_sliding_window(c_window)

    # Create list of positioned_loop_candidate objects
    loop_ensemble = []
    for ns in n_subframes:
        for cs in c_subframes:
            # Do not allow two single-residue frames--minimum three residues to align
            if len(ns + cs) == 2:
                continue

            # Create positioned_loop_candidate
            plc = positioned_loop_candidate(ns, cs, query_pose, subject_pose, matched_loop)
            
            if report:
                report.write(print_loop_candidate(plc))
            
            # Add to list
            loop_ensemble.append(plc)

    return loop_ensemble


def find_proximity_to_substrate(substrate_coords, loop_coords):
    """
    Finds all CA-CA distances between the coordinates of the substrate peptide
    and the candidate loop, if it were to be spliced into the query structure.
    Determines the cosest proximity and a list of residues in the peptide within 
    8A of the loop. 
    """
    # Finding close residues
    closest_distance = 1000 # Arbitrary large number
    nearby_substrate_residues = []
    for n, sr in enumerate(substrate_coords):
        for lr in loop_coords:
            substrate_dist = get_distance(sr, lr)
            if substrate_dist < closest_distance:
                closest_distance = substrate_dist
            if substrate_dist <= 8:
                if (n + 1) not in nearby_substrate_residues:
                    nearby_substrate_residues.append(n + 1)

    return closest_distance, nearby_substrate_residues


def CA_proximity_table(pose_1, selector_1, pose_2, selector_2):
    """
    Generates a CA distance table between selected residues that can be from 
    different poses.
    """
    # Get residue lists for table
    s1_res = selector_to_list(pose_1, selector_1)
    s2_res = selector_to_list(pose_2, selector_2)
    
    # Initialize table
    dist_table = np.zeros((len(s1_res), len(s2_res)))
    
    # Populate table
    for n1, i1 in enumerate(s1_res):
        ca1_coords = find_res_coords(pose_1, i1)
        for n2, i2 in enumerate(s2_res):
            ca2_coords = find_res_coords(pose_2, i2)
            dist_table[n1, n2] = get_distance(ca1_coords, ca2_coords)
            
    return dist_table


class positioned_loop_candidate():
    def __init__(self, n_frame, c_frame, query_pose, subject_pose, matched_loop):

        self.N_outside_overlap_residue = n_frame[0]
        self.N_splice_residue = n_frame[-1]
        self.C_splice_residue = c_frame[0]
        self.C_outside_overlap_residue = c_frame[-1]
        self.N_overlap_size = len(n_frame)
        self.C_overlap_size = len(c_frame)

        #
        q_list = get_residue_list(n_frame + c_frame, res_type='pose', target='query')
        s_list = get_residue_list(n_frame + c_frame, res_type='pose', target='subject')
        self.rmsd = align_protein_sections(query_pose, q_list, 
                Pose(subject_pose), s_list)

        # Identify clashes
        main_prot_select = ChainSelector('A')
        q_loop_n = get_ar_property(matched_loop.n_nearest_match, 'q', 'po')
        q_loop_c = get_ar_property(matched_loop.c_nearest_match, 'q', 'po')
        q_loop_select = ResidueIndexSelector('{}-{}'.format(q_loop_n, q_loop_c))
        template_select = AndResidueSelector()
        template_select.add_residue_selector(main_prot_select)
        template_select.add_residue_selector(NotResidueSelector(q_loop_select))
        s_loop_n = get_ar_property(matched_loop.n_nearest_match, 's', 'po') + 1
        s_loop_c = get_ar_property(matched_loop.c_nearest_match, 's', 'po') - 1
        s_loop_select = ResidueIndexSelector('{}-{}'.format(s_loop_n, s_loop_c))
        clash_table = CA_proximity_table(query_pose, template_select, 
            subject_pose, s_loop_select)
        self.clash_count = len(clash_table[clash_table < 4])

        # Identify substrate proximity
        subst_select = ChainSelector('B')
        substrate_table = CA_proximity_table(subject_pose, s_loop_select, 
            query_pose, subst_select)
        self.contact_count = len(substrate_table[substrate_table < 9])
        close_subst = np.where(substrate_table.min(0) <= 9)[0]
        self.close_substrate_residues = list(close_subst + 1) # +1 to 1-index

        # 
        self.is_near_target = bool(self.close_substrate_residues)
        self.is_low_clash = self.clash_count < 2


def pick_positioned_loop(candidates_list, selection='random', filter=True):
    """
    Pick a loop candidate out of a list of positioned_loop_candidate objects.
    Candidate can be selected on the basis of most substrate contacts ('contact'), 
    best RMSD to aligning residues at splice site ('rmsd'), or a list of substrate 
    residues that must be nearby (ex: [2,3,4,5,6]), or at random (default). If  
    multiple candidates are of equal quality, picks at random. By default, all  
    candidates with severe clashes or which are not near the substrate are filtered 
    out. If no acceptable loop is found, returns None.
    """
    # Validate options
    selection = selection.lower()
    if selection not in ['contact', 'rmsd', 'random']:
        if not isinstance(selection, list):
            print('Invalid selection')
            assert False

    # Filter
    if filter:
        viable_candidates = [c for c in candidates_list if c.is_near_target and c.is_low_clash]
        if len(viable_candidates) == 0:
            return None
    else: 
        viable_candidates = candidates_list

    # 
    if selection == 'random':
        return choice(viable_candidates)

    #
    if selection == 'contact':
        most_contacts = max([c.contact_count for c in viable_candidates])
        best_candidates = [c for c in viable_candidates if c.contact_count == most_contacts]

    if selection == 'rmsd':
        min_rmsd = min([c.rmsd for c in viable_candidates])
        best_candidates = [c for c in viable_candidates if c.rmsd == min_rmsd]

    if isinstance(selection, list):
        most_targets = max([len(list(set(c.contact_count) & set(selection))) for c in viable_candidates])
        if most_targets == 0:
            return None
        best_candidates = [c for c in viable_candidates if len(list(set(c.contact_count) & set(selection))) == most_targets]

    return choice(best_candidates)


def print_loop_candidate(positioned_loop_candidate, target='s', attribute='pdb'):
    """

    """
    # Template
    template = '{:<6}' * 9 + '{:<12}' + '{:<6}' * 2 + '\n'

    # Collect output data
    data_out = []
    for attr in ['N_outside_overlap_residue', 'N_splice_residue', 
        'C_splice_residue', 'C_outside_overlap_residue']:
        res = getattr(positioned_loop_candidate, attr)
        data_out.append(get_ar_property(res, target, attribute))

    for attr in ['N_overlap_size', 'C_overlap_size']:
        data_out.append(getattr(positioned_loop_candidate, attr))

    data_out.append(round(positioned_loop_candidate.rmsd,3))

    for attr in ['clash_count', 'contact_count']:
        data_out.append(getattr(positioned_loop_candidate, attr))

    close_res = positioned_loop_candidate.close_substrate_residues
    data_out.append(';'.join(close_res))

    for attr in ['is_near_target', 'is_low_clash']:
        data_out.append(getattr(positioned_loop_candidate, attr))

    return template.format(*data_out)


# Generating Dali files

# Running Dali on list
# 
# cd /mnt/c/Users/jhlub/Documents/dali_data/
# 
# ~/DaliLite.v5/bin/import.pl --pdbfile ../khare_lab/loop_swap_protocol/hcv/hcv.pdb --pdbid 0000 -dat DAT --clean
# 
# echo 0000A > hcv_query.list
# 
# ~/DaliLite.v5/bin/dali.pl --query hcv_query.list --db ../khare_lab/loop_swap_protocol/hcv/second_iteration_through_50_list.txt --dat1 DAT --dat2 DAT --clean 2> err
# 
# Output was 0000A.txt. Renamed to hcv_matches.txt and moved it to the hcv folder. Only 2hntC and 5mznA had Z-scores < 2.5.

# Generating HCV library

def retrieve_dali_pdb(pdb_id, pdb_dir):
    """
    For a given PDB ID, or PDB ID + Chain, determines the PDB path in my
    dali_data folder.
    """
    # If PDB ID includes chain, shrink to just PDB ID
    if len(pdb_id) > 4:
        pdb_id = pdb_id
        
    # Determine Dali folder for that PDB
    folder = pdb_id[1:3]
    
    # Generate file path
    return os.path.join(pdb_dir, '{}/pdb{}.ent.gz'.format(folder, pdb_id))


def unzip_gz(infile, outfile):
    """ Unzips a gzip file """
    import gzip
    import shutil

    with gzip.open(infile, 'rb') as f_in:
        with open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    return


def make_clean_pose(pdb_file, clean=True):
    """ 
    For a given PDB file, first checks if the file has a .gz extension
    and extracts it to an unzipped PDB file in the current directory. 
    Then runs the cleanATOM script on the PDB to strip out all non-ATOM lines, 
    and creates a .clean.pdb version in the local directory. Loads the PDB as 
    a pose, and then, by default, deletes the extracted and cleaned files 
    before returning the clean pose.
    """
    from pyrosetta.toolbox import cleanATOM

    bname = os.path.basename(pdb_file)
    
    # gzip extraction
    if pdb_file[-3:] == '.gz':
        extract = True
        unzipped_pdb = bname.replace('.gz', '.pdb')
        unzip_gz(pdb_file, unzipped_pdb)
    else:
        extract = False
        unzipped_pdb = pdb_file
    
    # Cleaning
    clean_pdb = os.path.basename(unzipped_pdb).replace('.pdb', '.clean.pdb')
    cleanATOM(unzipped_pdb, clean_pdb)
    
    # Load pose
    clean_pose = pose_from_pdb(clean_pdb)
    
    # Cleaning up temporary files
    if clean:
        os.remove(clean_pdb)
        if extract:
            os.remove(unzipped_pdb)
        
    return clean_pose


def get_single_chain_pose(pdb_id_chain, 
                          pdb_dir='/mnt/c/Users/jhlub/Documents/dali_data/pdb/'):
    """
    For a given 5-letter code corresponding to a 4-letter PDB ID 
    and a 1-letter chain ID, return a pose of just the appropriate chain.
    """
    assert isinstance(pdb_id_chain, str)
    assert len(pdb_id_chain) == 5
    
    pdb = pdb_id_chain[:4]
    chain = pdb_id_chain[-1].upper()
    
    pdb_path = retrieve_dali_pdb(pdb, pdb_dir)
    
    full_pose = make_clean_pose(pdb_path)
    
    for c in range(1, full_pose.num_chains() + 1):
        cstart = full_pose.chain_begin(c)
        cend = full_pose.chain_end(c)
        cid = full_pose.pdb_info().pose2pdb(cstart).split()[1]
        if cid == chain:
            chain_pose = Pose(full_pose, cstart, cend)
            break

    return chain_pose


def collect_list_of_homologs(dali_file):
    """
    Reads in a list of hits from a text file of Dali output 
    (can include alignments, etc.) and extracts just a list of PDB files 
    that have Z-scores of 2.5 or above.
    """
    with open(dali_file, 'r') as r:
        full_dali_summary = r.readlines()

    pdb_list = []
    for protein_neighbor_line in full_dali_summary[3:]: # Skip three header lines
        # Stop when getting to the alignment section
        if protein_neighbor_line == "\n":
            break

        # Identify PDB and chain, convert xxxx-C to xxxxC (C is chain)
        pnl_as_list = protein_neighbor_line.split()
        pdb = pnl_as_list[1].replace('-', '')
        
        # Collect only PDBs that have Z-scores >= 2.5
        zscore = float(pnl_as_list[2])
        if zscore >= 2.5:
            pdb_list.append(pdb)

    pdb_list.sort()
    
    return pdb_list


def collect_loop_library(query_name, query_pdb, pdb_list,
                         dali_hits_file, dali_alignment_file, 
                         catalytic_map, loop_map,
                         pdb_dir, subjects_format=0, 
                         pickle_file=None, dump_file=None):
    """
    Builds a loop library from a query PDB and Dali output, which can either
    be two separate files (list of hits and alignments) or a single file. 
    Also requires a directory of subject PDB files, assumed to be in the 
    format that Dali uses when syncing with the PDB. (Set subjects_format to 
    something else if you have already stripped the PDBs down to the chain 
    of interest.) Optionally writes a pickle file and data dump file. Returns 
    a list of any PDBs that failed to load or add properly to the library.
    """
    # Open query pose
    query_pose = pose_from_pdb(query_pdb)
    
    # Opening files
    if pickle_file:
        out_pickle = open(pickle_file, 'wb')
    if dump_file:
        out_dump = open(dump_file, 'w')
    else:
        out_dump = False
    
    # Create dictionary of subject PDBs and their loops
    homolog_loops_dict = {}
    fail_list = []
    for n, subj_name in enumerate(pdb_list):
        try:
            print('{:<5}/{:<5}:\t{:>10}'.format(n, len(pdb_list), subj_name))
            
            # Get subject pose, cleaning up from Dali if necessary
            if subjects_format == 0:
                subj_pose = get_single_chain_pose(subj_name, pdb_dir=pdb_dir)
            else:
                subj_pose = pose_from_pdb(join(pdb_dir, subj_name))
            
            # Map loops for the subject
            pinf = protease_info(query_name, subj_name, 
                                 dali_hits_file, dali_alignment_file, 
                                 query_pose, subj_pose, 
                                 catalytic_map, loop_map, 
                                 report=out_dump)
            homolog_loops_dict[subj_name] = pinf

            # Pickle results
            if pickle_file:
                out_dict = {subj_name: pinf}
                pickle.dump(out_dict, out_pickle)
                
        except:
            fail_list.append(subj_name)

    print('Failed:', len(fail_list))
    for i in fail_list:
        print('\t', i)

    # Closing files
    if pickle_file:
        out_pickle.close()
    if dump_file:
        out_dump.close()
    
    return homolog_loops_dict, fail_list


def unpickle_serial_dicts(pickle_file):
    """
    Takes a pickle file with a set of individually pickled dicts 
    (which is how collect_loop_library saves its data) and combines them
    into one dict.
    """
    serial_dict = {}
    
    with open(pickle_file, 'rb') as r:
        while True:
            try:
                unpickled_dict = pickle.load(r)
                serial_dict.update(unpickled_dict)
            except:
                break
    
    return serial_dict


def collect_valid_loops(homolog_loops_dict, 
                        summary_filename=None, pdb_dir=None):
    """
    A homolog_loops_dict includes many homologs, each with multiple loops.
    Extracts a list of loops that are candidates for insertion. If a 
    summary name is given, writes a summary of selected loops. Summary can 
    include a PyMOL command to load the loop if a directory is given for the
    source PDB files.
    """
    # Initialize collection list
    target_loops = []
    
    for prot in homolog_loops_dict.values():
        for lm in prot.loop_maps.values():
            if lm.is_possible_target:
                target_loops.append(lm)
                
    target_loops.sort(key=lambda x:(x.loop_name, x.loop_size_change, x.loop_source))
    
    if summary_filename:
        write_selected_loops_summary(target_loops, summary_filename, pdb_dir)
    
    return target_loops


def write_selected_loops_summary(loops_list, summary_filename, pdb_dir):
    """ 
    Writes a summary file for selected loops from collect_valid_loops.
    Output will be as a CSV, though the PyMOL command will be broken up 
    as a result. Summary can include a PyMOL command to load the loop if 
    a directory is given for the source PDB files. 
    
    Note: Command stil needs to be tweaked to accommodate CA-only alignmnt 
    and multiple aligned residues.
    """
    # Make header line
    loop_header = ['Loop', 'Source', 'Length', 'Length Change', 
                    'Distance to substrate', 'Close substrate residues',
                    'Query PDB Start', 'Query PDB End', 
                    'Subject PDB Start', 'Subject PDB End', 
                    'Query pose Start', 'Query pose End', 
                    'Subject pose Start', 'Subject pose End']
    if pdb_dir:
        loop_header.append('Pymol Command')

        # Make PyMOL command
        pymol_template =  "load {14}; "
        pymol_template += "create {0}_{1}, {1} and res {8}-{9}; "
        pymol_template += "delete {1}; "
        pymol_template += "pair_fit {0}_{1}///{8}+{9}/C+CA+N, tev///{6}+{7}/C+CA+N; "
        pymol_template += "hide everything, {0}_{1} and res {8}+{9}; "

    # Write summary file
    with open(summary_filename, 'w') as w:
        w.write(','.join(loop_header) + '\n')
        
        # Collect loop attributes for report
        for loop in loops_list:
            out_line = []
            out_line.append(loop.loop_name)                             # 0
            out_line.append(loop.loop_source)                           # 1
            out_line.append(loop.subject_loop_size)
            out_line.append(loop.loop_size_change)
            out_line.append(loop.closest_residue_distance)
            out_line.append(';'.join([str(i) for i in loop.close_substrate_residues]))
            out_line.append(loop.N_splice_residue.query_pdb_number)     # 6
            out_line.append(loop.C_splice_residue.query_pdb_number)     # 7
            out_line.append(loop.N_splice_residue.subject_pdb_number)   # 8
            out_line.append(loop.C_splice_residue.subject_pdb_number)   # 9
            out_line.append(loop.N_splice_residue.query_pose_number)
            out_line.append(loop.C_splice_residue.query_pose_number)
            out_line.append(loop.N_splice_residue.subject_pose_number)
            out_line.append(loop.C_splice_residue.subject_pose_number)
            out_line = [str(i) for i in out_line]

            if pdb_dir:
                pdb_path = retrieve_dali_pdb(loop.loop_source, pdb_dir) # 14
                out_line.append(pdb_path)
                com = com_temp.format(*out_line)
                out_line[14] = com
            
            w.write(','.join(out_line) + '\n')
            
    return


def collect_loops_subset_by_substrate_res(target_loops, substrate_res, 
    summary_filename=None, pdb_dir=None):
    """
    Takes a list of loop objects and a list or range of substrate residues
    and returns a subset of the loops list that includes only loops that 
    interact with those substrate residues.
    """
    # Collect list of nearby loops
    candidate_loops = []
    for loop in target_loops:
        for csr in loop.close_substrate_residues:
            if csr in substrate_res:
                candidate_loops.append(loop)
                break
    
    if summary_filename:
        write_selected_loops_summary(candidate_loops, summary_filename, pdb_dir)
        
    return candidate_loops


# Grafting in loops

def pick_random_loop(candidate_loops):
    """ 
    From a list of candidate loops, selects one loop at random that
    is a possible target.
    """    
    selected_loop = None
    while not selected_loop:
        # Randomly select a loop
        candidate = choice(candidate_loops)
        
        # Check that loop is a suitable target
        if candidate.is_possible_target:
            selected_loop = candidate
            
        else:
            candidate_loops.remove(candidate)
            
    return selected_loop


def extract_loop_subpose(pose, loop, dump_loop_pdb=None,
                        template_pose=None, aligned_residues=None):
    """
    Given a pose and matched_loop object, extract a subpose based on the 
    splice sites identified in the matched_loop. Can dump the loop pose 
    optionally--if option is True, dumps in the working directory. If a 
    string is given, that will be the output directory. If an align_template 
    pose and aligned_residues list are provided, will align poses before 
    doing the extraction.
    """
    # Pre-align
    if template_pose and aligned_residues:
        prealign_by_dali(template_pose, pose, aligned_residues)
    
    # Determine terminal residues
    extract_start = get_ar_property(loop.N_outside_overlap_residue, 's', 'po')
    extract_end = get_ar_property(loop.C_outside_overlap_residue, 's', 'po')
    
    # Make extract pose
    loop_pose = Pose(pose, extract_start, extract_end)
    
    # Dumping
    if dump_loop_pdb:
        extract_name = 'loop{}_{}.pdb'.format(loop.loop_name, 
                                       loop.loop_source)
        if isinstance(dump_loop_pdb, str):
            loop_pose.dump_pdb(join(dump_loop_pdb, extract_name))
        else:
            loop_pose.dump_pdb(extract_name)
            
    return loop_pose


def graft_candidate_loop(template_pose, loop_pose, loop):
    """
    Creates a pose based on a template, but with a loop swapped in
    based on a matched_loop object
    """
    # Determine insert sites
    n_splice = get_ar_property(loop.N_splice_residue, 'q', 'po')
    c_splice = get_ar_property(loop.C_splice_residue, 'q', 'po')
    
    # Determine overlaps
    n_overlap = loop.N_overlap_size
    c_overlap = loop.C_overlap_size
    
    # Set up CCDEndsGraftMover
    ccdgm = CCDEndsGraftMover()
    ccdgm.set_insert_region(n_splice, c_splice)
    ccdgm.set_piece(loop_pose, n_overlap, c_overlap)
    
    # Applying
    swapped_pose = Pose(template_pose)
    ccdgm.apply(swapped_pose)
    
    return swapped_pose


# Running after 2/26/20 rebuild

def harsh_collect_loop_library(query_name, query_pdb, pdb_list,
                         dali_hits_file, dali_alignment_file, 
                         catalytic_map, loop_map,
                         pdb_dir, subjects_format=0, 
                         pickle_file=None, dump_file=None):
    # Open query pose
    query_pose = pose_from_pdb(query_pdb)
    
    # Opening files
    if pickle_file:
        out_pickle = open(pickle_file, 'wb')
    if dump_file:
        out_dump = open(dump_file, 'w')
    else:
        out_dump = False
    
    # Create dictionary of subject PDBs and their loops
    homolog_loops_dict = {}
    for n, subj_name in enumerate(pdb_list):
        print('{:<5}/{:<5}:\t{:>10}'.format(n, len(pdb_list), subj_name))

        # Get subject pose, cleaning up from Dali if necessary
        if subjects_format == 0:
            subj_pose = get_single_chain_pose(subj_name, pdb_dir=pdb_dir)
        else:
            subj_pose = pose_from_pdb(join(pdb_dir, subj_name))

        # Map loops for the subject
        pinf = protease_info(query_name, subj_name, 
                             dali_hits_file, dali_alignment_file, 
                             query_pose, subj_pose, 
                             catalytic_map, loop_map, 
                             report=out_dump)
        homolog_loops_dict[subj_name] = pinf

        # Pickle results
        if pickle_file:
            out_dict = {subj_name: pinf}
            pickle.dump(out_dict, out_pickle)

    # Closing files
    if pickle_file:
        out_pickle.close()
    if dump_file:
        out_dump.close()
        
    return homolog_loops_dict


def parse_args():
    info=''
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument("-qn", "--query_name", default='HCV', 
        type=str, help="ID of template protease")
    parser.add_argument("-qpdb", "--query_pdb", default='hcv/hcv.pdb', 
        type=str, help="PDB file of template protease")
    parser.add_argument("-pdb_dir", "--subject_pdb_dir", default='/mnt/c/Users/jhlub/Documents/dali_data/pdb', 
        type=str, help="Directory location for subject PDB files")
    parser.add_argument("-spdb", "--subject_pdbs", default=None, 
        help="Input subject PDBs. By default, generates a list from the Dali \
        hits file. If a text file is given where each line is a PDB, reads that.")
    parser.add_argument("-dh", "--dali_hits", default='hcv/hcv_matches.txt', 
        type=str, help="Dali homolog list")
    parser.add_argument("-da", "--dali_alignments", default='hcv/hcv_matches.txt', 
        type=str, help="Dali alignments list")
    parser.add_argument("-pckl", "--pickle_file", default='loop_database.pkl',
        help="To what pickle file should the database write? Strongly recommended to use this feature.")
    parser.add_argument("-dd", "--data_dump", 
        help="To what text file should the database generation output be written?")
    parser.add_argument("-f", "--parallel_fraction", nargs=2, type=int,
        help="For parallelization, put in two values--the first is this job's place in the group, \
        the second is the total number of processors. So 2 5 would take the second 20pct of the total.")
    
    return parser.parse_args()


def main(args):
    print(args)
    init('-mute all -ignore_zero_occupancy false')

    # Global data
    hcv_cat_res = {'H':72, 'A':96, 'N':154}
    hcv_map     = {'N':  range(  1,   4),
                    1:   range( 11,  21),
                    2:   range( 25,  47),
                    3:   range( 53,  57),
                    4:   range( 64,  66),
                    5:   range( 71,  89),
                    6:   range( 93,  96),
                    7:   range(102, 118),
                    8:   range(124, 126),
                    9:   range(136, 138),
                    10:  range(142, 155),
                    11:  range(161, 163),
                    12:  range(176, 178),
                   'C':  range(187, 197)} 
    
    if args.subject_pdbs:
        with open(args.subject_pdbs, 'r') as r:
            subj_pdbs = r.readlines()
            pdb_list = [p.strip() for p in subj_pdbs]
    else:
        pdb_list = collect_list_of_homologs(args.dali_hits)
        print(len(pdb_list))

    
    if args.parallel_fraction:
        fraction_len = round(len(pdb_list)/args.parallel_fraction[1])
        fraction_begin = fraction_len * (args.parallel_fraction[0] - 1)
        if args.parallel_fraction[0] == args.parallel_fraction[1]:
            fraction_end = len(pdb_list) + 1
        else:
            fraction_end = fraction_len * (args.parallel_fraction[0])
        pdb_list = pdb_list[fraction_begin: fraction_end]
    print(len(pdb_list))

    db, fail_list = collect_loop_library(args.query_name, args.query_pdb, pdb_list,
                         args.dali_hits, args.dali_alignments, 
                         hcv_cat_res, hcv_map,
                         args.subject_pdb_dir, subjects_format=0, 
                         pickle_file=args.pickle_file)
    print(fail_list)


    

if __name__ == '__main__':
    args = parse_args()
    main(args)

