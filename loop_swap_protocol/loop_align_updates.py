#!/usr/bin/python
from Bio import pairwise2
from glob import glob
from math import sqrt
import numpy as np
import pickle
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import rmsd_atoms
from pyrosetta.rosetta.core.scoring import superimpose_pose
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import \
    PerResidueRMSDMetric

################################################################################
# General utility functions

def get_distance(c1, c2):
    """ Returns the distance between two Rosetts XYZ coordinate vectors"""
    dist = sqrt((c2.x - c1.x) ** 2 + (c2.y - c1.y) ** 2 + (c2.z - c1.z) ** 2)
    return dist


def find_res_coords(pose, resnum, atom_type='CA'):
    """ For a given pose and residue number, returns the coordinates of a 
    specified atom type """
    residue = pose.residue(resnum)
    res = residue.atom(atom_type)
    return res.xyz()


def list_pose_coords(pose, atom_type='CA'):
    """ For a given pose, list all CA coordinates """
    # Initialize list of CA coordinates
    ca_coords = []

    # Populating list of CA coordinates 
    for res in range(1, pose.total_residue() + 1):
        res_ca = find_res_coords(pose, res, atom_type=atom_type)
        ca_coords.append(res_ca)

    return ca_coords


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


def align_protein_sections(pose_1, selector_1, pose_2, selector_2, mode='CA'):
    """
    Aligns selected regions of two poses, superimposing the second pose onto 
    the first, based on RMSD. Returns the RMSD value. Input selections can
    either be Rosetta selector objects or lists of residue by pose number, 
    which will be made into a selector. By default, uses CA RMSD. Can also 
    use full backbone.
    """
    # Verify mode is acceptable
    assert mode in ['CA', 'BB']
    
    # Set RMSD type based on mode
    if mode == 'CA':
        rmsd_type = rmsd_atoms.rmsd_protein_bb_ca
    if mode == 'BB':
        rmsd_type = rmsd_atoms.rmsd_protein_bb_heavy
    
    #print(selector_1, selector_2)
    # If a list is given, make a selector
    if isinstance(selector_1, list):
        select_str_1 = ','.join([str(i) for i in selector_1])
        selector_1 = ResidueIndexSelector(select_str_1)

    if isinstance(selector_2, list):
        select_str_2 = ','.join([str(i) for i in selector_2])
        selector_2 = ResidueIndexSelector(select_str_2) 

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
    return round(total_b / 3, 3)


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


#def find_max_distance(coords_list):
    #   """
    #   Given a list of coordinates, determines the maximum distance between any 
    #   pair.
    #   """
    #   max_distance = 0
    #
    #   # Comparing all coordinate pairs
    #   for n, coord in enumerate(coords_list):
    #       for partner in coords_list[n + 1:]:
    #           distance = get_distance(coord, partner)
    #
    #           if distance > max_distance:
    #               max_distance = distance
    #
    #   return max_distance

################################################################################
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
        query_pose, subject_pose, catalytic_residues, structure_map):

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
        self.loop_maps = {}

        self.auto_calculate(dali_file, align_file, query_pose, subject_pose, 
            catalytic_residues, structure_map)

    def auto_calculate(self, dali_file, align_file, query_pose, subject_pose, 
        catalytic_residues, structure_map):
        """ Runs all update methods """
        if dali_file:
            self.update_dali_info(dali_file)

        if all([align_file, query_pose, subject_pose]):
            self.update_aligned_residues(align_file, query_pose, subject_pose)

        if all([self.aligned_residues, catalytic_residues]):
            self.update_catalytic_residues(catalytic_residues)

        if all([self.aligned_residues, query_pose, subject_pose, 
            structure_map]):
            self.update_loop_maps(structure_map, query_pose, subject_pose)

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

    def update_aligned_residues(self, align_file, query_pose, subject_pose):
        """
        Updates aligned_residues property with a list populated from a Dali 
        alignmant file and query and subject poses. The latter are necessary 
        to determine pose numbers and B-factors.
        """
        self.aligned_residues = map_aligned_residues(align_file, 
            self.subject_name, query_pose, subject_pose)
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

    def update_loop_maps(self, structure_map, query_pose, subject_pose):
        """
        Updates loop_maps attribute with a list of matched_loop objects, based 
        on a structure map and query and subject poses. The structure map should 
        be a dict of the form {'N': range(<N-term_start>, <N-term_end>), 
        1:range(<loop_1_satrt>,<loop_1_end>), ..., 'C': range(<C-term_start>, 
        <C-term_end>)}. The ranges are assumed to be in PDB numbers.
        """
        self.loop_maps = map_structure_elements(structure_map, query_pose, 
            subject_pose, self.aligned_residues, self.subject_name)
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
        dali_info['Z_score'] = summary_line[1]
        dali_info['rmsd'] = summary_line[2]
        dali_info['lali'] = summary_line[3]
        dali_info['nres'] = summary_line[4]
        dali_info['pID'] = summary_line[5]
        dali_info['description'] = ' '.join(summary_line[6:])

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
                    end_block = n 
                    break
        
                # Find beginning of block, where start line includes name
                if subject in sa.upper():
                    begin_block = n

        # Extracting relevant text block
        alignment_block = seq_aligns[begin_block:end_block]

        return alignment_block


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
                end_block = n 
                break
    
            # Find beginning of block, where start line includes name
            if subject in sa.upper():
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


def get_posenums_from_dali_str(pose, dali_string):
    """
    For a given pose and matching string from a Dali alignment, returns a list 
    of pose numbers. Where there ar gaps, the list will include a None entry. 
    Uses Biopython's global alignment with maximum similarity function, with 
    scores from https://towardsdatascience.com 
    Matching characters: +2
    Mismatching character: -1
    Opening a gap: -0.5
    Extending a gap: -0.1
    Assumes that if the pose has multiple chains, the first is the chain that 
    is aligned. 
    """ 
    # If pose has multiple chains, only take first one
    pose_chain_1 = pose.split_by_chain()[1]

    # Get sequence strings from pose and Dali
    ps = pose_chain_1.sequence().upper()    # Pose sequence
    ds = dali_string.upper()        # Dali sequence

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

    return pnlist[:len(ds)] 


def map_aligned_residues(align_file, subject, query_pose, subject_pose):
    """
    Feed alignment data into a list of aligned_residue objects, each with 
    corresponding information about the position in both query and the subject 
    protease being analyzed.
    """
    # Get alignment 
    alignment = get_dali_alignment(align_file, subject)

    # Get pose number lists, add to alignment
    alignment['query_pose_numbers'] = get_posenums_from_dali_str(query_pose, 
        alignment['query_sequence'])
    alignment['subject_pose_numbers'] = get_posenums_from_dali_str(subject_pose, 
        alignment['subject_sequence'])

    # Initialize aligned residues list
    aligned_residues = []

    # Loop through each residue in the alignment, adding aligned_residue
    # objects to the list for each
    for i in range(alignment['length']):
        # Check that at least one pose-numbered residue is present
        q_res = alignment['query_pose_numbers'][i]
        s_res =  alignment['subject_pose_numbers'][i]
        if not any([q_res, s_res]):
            continue

        # Populating aligned_residue object
        a_residue = aligned_residue(i, alignment, query_pose, subject_pose)

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
    aligned_residues, subject_name, mode='pdb', target='query'):
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
    assert mode.lower() in ['pdb', 'pose']

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
            n_bound = None
            c_bound = structure_map[1][0] - 1 
            # ^ = Everything up to first res of first loop
        elif loop == 'C': # Edge case for C-terminal region (not a loop)
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
            n_bound, c_bound, mode=mode, target=target)[1]

        # Make matched_loop object and add it to the dict
        loop_map = matched_loop(subject_name, loop, ar_subset, 
            structure_map[loop], query_pose, subject_pose)
        loop_maps[loop] = loop_map

    return loop_maps

################################################################################
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
    def __init__(self, index, alignment, query_pose, subject_pose):
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

        self.auto_populate(index, alignment, query_pose, subject_pose)

    def auto_populate(self, index, alignment, query_pose, subject_pose):
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


def partition_aligned_residues_list(aligned_residues, n_bound, c_bound, 
    mode='pdb', target='query'):
    """
    Returns a subsets of a list of aligned_residue objects, partitioned with a  
    given set of boundaries. The boundaries can be either PDB or pose numbers.  
    The bounds can correspond to residues in either the query protein or the 
    subject. Returns a list of three lists: the set before, the set between, 
    and the set after. The bounds are included in the middle set.
    """
    # Make sure mode is either pdb or pose
    assert mode.lower() in ['pdb', 'pose']

    # Make sure target is either query or subject
    assert target.lower() in ['query', 'subject']

    # Loop through list to find boundary residues for given mode and target
    for n, ar in enumerate(aligned_residues):
        if mode.lower() == 'pdb':
            if target.lower() == 'query':
                if ar.query_pdb_number == n_bound:
                    start_index = n
                if ar.query_pdb_number == c_bound:
                    end_index = n
            if target.lower() == 'subject':
                if ar.subject_pdb_number == n_bound:
                    start_index = n
                if ar.subject_pdb_number == c_bound:
                    end_index = n
        if mode.lower() == 'pose':
            if target.lower() == 'query':
                if ar.query_pose_number == n_bound:
                    start_index = n
                if ar.query_pose_number == c_bound:
                    end_index = n
            if target.lower() == 'subject':
                if ar.subject_pose_number == n_bound:
                    start_index = n
                if ar.subject_pose_number == c_bound:
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
    # Initialize collection dict
    matches = {'nearest_match': None, 'nearest_b_match': None,
               'farthest_match': None, 'farthest_b_match': None}

    # Scan through list of aligned residues
    for ar in aligned_residues:
        # Ignore all residues that do not align until first match found
        if ar.residues_align:

            # Continue updating farthest_match through list
            matches['farthest_match'] = ar

            # Only capture first match as nearest_match
            if not matches['nearest_match']:
                matches['nearest_match'] = ar

            # Similar process, but limited to B-sheet residues
            if ar.subject_sec_struct == 'E':
                matches['farthest_b_match'] = ar
                if not matches['nearest_b_match']:
                    matches['nearest_b_match'] = ar

        # After finding matches, stop searching at the first unmatched residue
        else:
            if matches['nearest_match']:
                break

    return matches


def get_residue_list(aligned_residues, mode='pdb', target='query'):
    """
    For a given list of aligned_residue objects, returns a list of residues of 
    the specified type. Default behavior is to return the PDB numbers of the 
    query protein, though the subject protein can be searched instead, and the 
    pose numbers or B factors can be searched instead.
    """
    # Make sure mode is either pdb or pose
    assert mode.lower() in ['pdb', 'pose', 'b']

    # Make sure target is either query or subject
    assert target.lower() in ['query', 'subject']

    # Initialize output list
    out_res_list = []

    # Loop through list and collect the desired attribute
    for ar in aligned_residues:
        if mode == 'pdb' and target == 'query':
            out_res_list.append(ar.query_pdb_number)

        if mode == 'pdb' and target == 'subject':
            out_res_list.append(ar.subject_pdb_number)

        if mode == 'pose' and target == 'query':
            out_res_list.append(ar.query_pose_number)

        if mode == 'pose' and target == 'subject':
            out_res_list.append(ar.subject_pose_number)

        if mode == 'b' and target == 'query':
            out_res_list.append(ar.query_b_factor)

        if mode == 'b' and target == 'subject':
            out_res_list.append(ar.subject_b_factor)

    return out_res_list

################################################################################
# Functions used by matched_loop class

class matched_loop():
    """
    Data storage structure for loops. When taking in a loop of a query 
    structure, finds the edges bordering it (usually B-sheets) and looks for 
    residue matches within the given boundaries, which should be the starts of
    the next loops. Input residues should use PDB (not pose) numbers.
    """
    def __init__(self, source, l_name, local_residues, l_range, query_pose, 
        subject_pose):

        # Loop source and name
        self.loop_source = source
        self.loop_name = l_name

        # Flanking residues -- update_matches
        self.n_nearest_match = None
        self.n_farthest_match = None
        self.c_nearest_match = None
        self.c_farthest_match = None
        self.flanking_matches_found = False

        # Best matched residues for loop swap -- update_splice_sites
            # Residues listed in sequence order N to C
        self.N_outside_overlap_residue = None
        self.N_splice_residue = None
        self.N_loop_end_residue = None
        self.C_loop_end_residue = None
        self.C_splice_residue = None
        self.C_outside_overlap_residue = None
            # Numerical properties of loop
        self.N_overlap_size = None
        self.C_overlap_size = None
        self.query_loop_size = None
        self.subject_loop_size = None
        self.loop_size_change = None

        # Loop RMSD is query and subject are the same size -- update_loop_rmsd
        self.rmsd = None

        # Average B factor for the loop -- update_b_factor
        self.b_factor = None

        # Check for discontinuities -- update_continuity
        self.is_continuous = None

        # Loop proximity to peptide substrate -- update_proximity_to_substrate
        self.closest_residue_distance = None
        self.close_substrate_residues = []

        # Evaluate whether loop is a suitable target -- update_suitability
        self.is_n_match = None
        self.is_c_match = None
        self.is_near_target = None
        self.is_not_domain = None
        self.is_different_from_original = None
        self.is_ordered = None
        self.is_possible_target = None

        self.auto_calculate(local_residues, l_range, query_pose, subject_pose)

    def auto_calculate(self, local_residues, l_range, query_pose, subject_pose):
        """ Runs all update methods """
        if all([local_residues, l_range]):
            res_sets = self.update_matches(local_residues, 
                l_range[0], l_range[-1])

        if all([res_sets, query_pose, subject_pose]):
            n_splice_range, c_splice_range = self.update_splice_sites(res_sets, 
            query_pose, subject_pose)
            self.update_loop_rmsd(query_pose, subject_pose)
            self.update_b_factor(local_residues)
            self.update_continuity(subject_pose)

        if all([query_pose.num_chains() == 2, query_pose, subject_pose, 
            n_splice_range, c_splice_range]):
            self.update_proximity_to_substrate(query_pose, subject_pose, 
                n_splice_range + c_splice_range)

        self.update_suitability()

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
            loop_start, loop_end, mode='pdb', target='query')

        # Collect matching residues, reversing the n-side set so as to count 
        # away from the loop.
        n_matches = find_nearest_matches(res_sets[0][::-1])
        c_matches = find_nearest_matches(res_sets[2])

        # Update properties
        self.n_nearest_match = n_matches['nearest_match']
        self.n_farthest_match = n_matches['farthest_match']
        self.c_nearest_match = c_matches['nearest_match']
        self.c_farthest_match = c_matches['farthest_match']

        # Indicate if matches were found
        if self.n_nearest_match and self.c_nearest_match:
            self.flanking_matches_found = True

        return res_sets
        
    def update_splice_sites(self, res_sets, query_pose, subject_pose):
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
            query_pose, subject_pose)
        
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
        self.query_loop_size = self.C_splice_residue.query_pose_number - \
            self.N_splice_residue.query_pose_number - 1
        self.subject_loop_size = self.C_splice_residue.subject_pose_number - \
             self.N_splice_residue.subject_pose_number - 1
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
        query_list = list(range(self.N_splice_residue.query_pose_number,
            self.C_splice_residue.query_pose_number + 1))
        subject_list = list(range(self.N_splice_residue.subject_pose_number,
            self.C_splice_residue.subject_pose_number + 1))

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
        first_loop_res = local_residues.index(self.N_splice_residue) + 1
        last_loop_res = local_residues.index(self.C_splice_residue)
        loop_residues = local_residues[first_loop_res:last_loop_res]

        # Collect list of B factors 
        b_factors = get_residue_list(loop_residues, mode='b', target='subject')

        # Exclude empty values
        b_factors_clean = [i for i in b_factors if isinstance(i, float)]

        # Take average backbone B factor for the loop
        b_average = np.average(b_factors_clean)

        # Update attribute
        self.b_factor = b_average

        return

    def update_continuity(self, subject_pose):
        """
        Checks whether the loop includes chain breaks. Updates the is_continuous 
        attribute.
        """
        # Nothing to update if no matching residues were identified
        if not self.flanking_matches_found:
            return

        # Create loop subpose
        subpose_start = self.N_outside_overlap_residue.subject_pose_number
        subpose_end = self.C_outside_overlap_residue.subject_pose_number
        subpose = Pose(subject_pose, subpose_start, subpose_end)

        # Check for breaks
        continuous, c_n_distances, break_sites = check_pose_continuity(subpose)

        # Announce identified breaks
        if not continuous:
            print('Chain broken:')
            for b in breaks:
                print('\t', subpose.pdb_info().pose2pdb(b).split()[0])

        # Update attribute
        self.is_continuous = continuous

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
        # Getting coordinates
        loop_start = self.N_splice_residue.subject_pose_number + 1
        loop_end = self.C_splice_residue.subject_pose_number - 1
        substrate_coords, loop_coords = get_substrate_and_loop_coords(
            query_pose, subject_pose, overlaps, loop_start, loop_end)

        # Determine closest proximity and list of nearby substrate residues
        closest_distance, nearby_substrate_residues = \
            find_proximity_to_substrate (substrate_coords, loop_coords)

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
            n_match_check = bool(self.N_splice_residue)

        # Check that there is a matched residue on the C-terminal side
        if self.loop_name == 'C':
            c_match_check = True
        else:
            c_match_check = bool(self.C_splice_residue)

        # Check that there are residues within range of substrate
        proximity_check = bool(self.close_substrate_residues)

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

        # Updating attributes
        self.is_n_match = n_match_check
        self.is_c_match = c_match_check
        self.is_near_target = proximity_check
        self.is_not_domain = res_count_check
        self.is_different_from_original = similarity_check
        self.is_ordered = b_factor_check
        self.is_possible_target = all([n_match_check, c_match_check, 
            proximity_check, res_count_check, similarity_check, b_factor_check, 
            self.is_continuous])

        return


def find_splice_sites(n_window, c_window, query_pose, subject_pose):
    """
    For a given pair of residue windows on the N-side and C-side of a loop, 
    break each window into subframes of length varied between 2 and the full 
    window size, unless the starting window is only of lenght 1, in which case 
    it is used directly. Then take those subframes from each and align based on  
    the residues in those subframes, and determine the best backbone RMSD score.
    The combination of residues that yields the lowest RMSD will be returned.
    """
    # List subframes within overlap windows
    if len(n_window) > 1:
        n_subframes = variable_sliding_window(n_window, min_size=2)
    else:
        n_subframes = variable_sliding_window(n_window)

    if len(c_window) > 1:
        c_subframes = variable_sliding_window(c_window, min_size=2)
    else:
        c_subframes = variable_sliding_window(c_window)

    # Initialize collection objects for best alignment
    best_rmsd = 1000 # Arbitrarily high
    best_n_set = n_subframes[-1] # Defaults to the full alignment region
    best_c_set = c_subframes[-1] # Defaults to the full alignment region

    # Loop through all subframes and align based on those residues, taking 
    # the set with the lowest RMSD
    for ns in n_subframes:
        for cs in c_subframes:
            # Get pose residues list for the residues to take RMSD
            q_list = get_residue_list(ns + cs, mode='pose', target='query')
            s_list = get_residue_list(ns + cs, mode='pose', target='subject')

            # Align these residues and get RMSD
            rmsd = align_protein_sections(query_pose, q_list, 
                Pose(subject_pose), s_list)

            # Compare RMSD to best, update best if lower
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_n_set = ns
                best_c_set = cs

    return best_n_set, best_c_set


class positioned_loop_candidate():
    def __init__(self, n_frame, c_frame, query_pose, subject_pose, 
        structure_map):

        self.rmsd = None
        self.ccd_graft_str = None


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

################################################################################



def main():
    # Initialize Rosetta
    init('-mute all, -ignore_zero_occupancy false')

    # Global data
    hcv_cat_res = {'H':72, 'A':96, 'N':154}
    tev_cat_res = {'H': 46, 'A': 81, 'N': 151} # In pose: H39, D74, C144

    # Map lists all loop regions, catalytic triad, all in PDB (not pose) numbers
    # Had to generate manually since loops might contain beta-fingers
    tev_map = {    'N':  range(  8,  18),
                    1:   range( 26,  28),
                    2:   range( 38,  40),
                    3:   range( 44,  54),
                    4:   range( 59,  63),
                    5:   range( 67,  74),
                    6:   range( 78,  82),
                    7:   range( 87, 109),
                    8:   range(117, 121),
                    9:   range(126, 139),
                    10:  range(143, 152),
                    11:  range(158, 161),
                    12:  range(172, 176),
                   'C':  range(182, 221)}

    htra1_map = {  'N':  range(160, 183),
                    1:   range(192, 198),
                    2:   range(209, 214),
                    3:   range(218, 226),
                    4:   range(233, 235),
                    5:   range(240, 241),
                    6:   range(247, 250),
                    7:   range(256, 276),
                    8:   range(284, 290),
                    9:   range(300, 317),
                    10:  range(320, 329),
                    11:  range(335, 337),
                    12:  range(349, 351),
                   'C':  range(357, 370)}

    hcv_map =   {  'N':  range(  1,   4),
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

    dali_file='aligned_pdbs/0000_dali_pdb90_tev.txt'
    align_file='aligned_pdbs/0000_seq_align.txt'

    pdb_list = glob('aligned_pdbs/*.pdb')
    pdb_list.sort()
    pdb_list.remove('aligned_pdbs/0000_master_pdb.pdb')

    db_dict = {}
    fail_list = []
    with open('alignment_data_dump.txt', 'w') as w:
        for i in pdb_list:
            print(i)
            subj_name = i.replace('aligned_pdbs/','').replace('.pdb', '')
            subj_pose = pose_from_pdb(i)
            try:
                pinf = protease_info(query_name, subj_name, dali_file, align_file, 
                                     query_pose, subj_pose, tev_cat_res, tev_map)
                db_dict[subj_name] = pinf
            except:
                fail_list.append(i)
                print(i, 'fail')

    print(len(fail_list))
    for i in fail_list:
        print('\t', i)

    outfile = 'protease_db.pkl'
    with open(outfile, 'wb') as o:
        pickle.dump(db_dict, o)


if __name__ == '__main__':
    main()

"""
pdb_list = glob('aligned_pdbs/*.pdb')
pdb_list.sort()
pdb_list.remove('aligned_pdbs/0000_master_pdb.pdb')
db_collect = []
for i in pdb_list:
    print(i)
    try:
        subj_name = i.replace('aligned_pdbs/','').replace('.pdb', '')
        subj_pose = pose_from_pdb(i)
        pinf = protein_loop_alignment.protease_info('TEV', tev, subj_name, subj_pose)
        pinf.auto_calculate()
        db_collect.append(pinf)
    except:
        print(i, 'failed')

outfile = 'protease_db.pkl'
with open(outfile, 'wb') as o:
    pickle.dump(db_collect, o)

with open('protease_database.csv', 'w') as w:
    header = ['query', 'subject', 'Z_score', 'rmsd', 'lali', 'nres', 'pID']
    header += ['cat_nucleophile', 'cat_his', 'cat_acid']
    for k in tev_map.keys():
        header += ['loop', 'length', 'potential_target', 'query_range', 'subject_range', 'reasons_rejected']
    w.write(', '.join(header) + '\n')

with open('protease_database.csv', 'a') as w:
    for i in db_collect:
        line_info = []
        line_info = [i.query_name, i.subject_name, i.Z_score, i.rmsd, i.lali, i.nres, i.pID]
        if i.nucleophile_res:
            line_info.append(i.nucleophile_type + str(i.nucleophile_res))
        else:
            line_info.append('None')
        if i.catalytic_his:
            line_info.append(i.catalytic_his_type + str(i.catalytic_his))
        else:
            line_info.append('None')
        if i.catalytic_acid:
            line_info.append(i.catalytic_acid_type + str(i.catalytic_acid))
        else:
            line_info.append('None')
        for k in i.loop_maps.keys():
            header += ['loop', 'length', 'swap_target', 'query_range', 'subject_range']
            line_info += [k, i.loop_maps[k].residue_count, i.loop_maps[k].is_possible_target]
            if i.loop_maps[k].is_possible_target:
                line_info.append('-'.join([str(x) for x in [i.loop_maps[k].query_N_splice_res, i.loop_maps[k].query_C_splice_res]]))
                line_info.append('-'.join([str(x) for x in [i.loop_maps[k].subject_N_splice_res, i.loop_maps[k].subject_C_splice_res]]))
                line_info += ['']
            else:
                line_info += ['', '']
                reject_reasons = []
                if not i.loop_maps[k].is_near_target:
                    reject_reasons.append('Distance')
                if not i.loop_maps[k].is_not_domain:
                    reject_reasons.append('Size')
                if not i.loop_maps[k].is_n_match:
                    reject_reasons.append('No N match')
                if not i.loop_maps[k].is_c_match:
                    reject_reasons.append('No C match')
                line_info.append('; '.join(reject_reasons))
        line_info = [str(i) for i in line_info]
        w.write(', '.join(line_info) + '\n')

"""

"""

only allow D&E for acids
"""

# Need to bring back terminal regions -- presently precluded by flanking_matches_found check
    # Check if each side is found, then have process look for adjacent loops that can add up to the whole
    # Every other? (sheets going opposite directions)
# Turn nobs on sliding window overlaps -- only picking single residues
# Warning: Pose includes residues beyond Dali alignment.
# Address loops that displace multiple loops. Ex: 1CU1 connects from start of loop 3 to end of loop 5
# Initial fail list: 23
#      aligned_pdbs/1CU1.pdb
#      aligned_pdbs/1EZX.pdb
#      aligned_pdbs/1FDP.pdb
#      aligned_pdbs/1FY1.pdb
#      aligned_pdbs/1KY9.pdb -- Worked before
#      aligned_pdbs/1PYT.pdb
#      aligned_pdbs/1RGQ.pdb
#      aligned_pdbs/1SGF.pdb
#      aligned_pdbs/1SHY.pdb
#      aligned_pdbs/1SVP.pdb
#      aligned_pdbs/1lvo.pdb
#      aligned_pdbs/2WV9.pdb
#      aligned_pdbs/2m9p.pdb -- Worked before
#      aligned_pdbs/3LKW.pdb
#      aligned_pdbs/3U1I.pdb
#      aligned_pdbs/4M9M.pdb
#      aligned_pdbs/4R8T.pdb
#      aligned_pdbs/5ILB.pdb
#      aligned_pdbs/5LC0.pdb
#      aligned_pdbs/5MZ4.pdb
#      aligned_pdbs/5WDX.pdb
#      aligned_pdbs/5t1v.pdb
#      aligned_pdbs/6BQJ.pdb
# Report is misplacing the following in verbose mode:
#   Residue alignments:
#                   Query                                   Subject                                 
#   Equal   Aligned AA_type PDB_no  Pose_no DSSP    Bfactor AA_type PDB_no  Pose_no DSSP    Bfactor 
# Make checks on all bad loops by failure type