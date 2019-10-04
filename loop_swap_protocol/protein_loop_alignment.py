#!/usr/bin/python
from pyrosetta import *
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.protocols.simple_moves import SuperimposeMover
from math import sqrt


init('-mute all')
#hcv_pose = pose_from_pdb('a_to_s_ly104_WT.pdb')
#hcv_cat_res = {72:'H', 96:'D', 154:'S'}
tev = pose_from_pdb('tev.pdb')
tev_seq = 'GESLFKGPRDYNPISSTICHLTNESDGHTTSLYGIGFGPFIITNKHLFRRNNGTLLVQSLHGVFKVKNTTTLQQHLIDGRDMIIIRMPKDFPPFPQKLKFREPQREERICLVTTNFQTKSMSSMVSDTSCTFPSSDGIFWKHWIQTKDGQCGSPLVSTRDGFIVGIHSASNFTNTNNYFTSVPKNFMELLTNQEAQQWVSGWRLNADSVLWGGHKVFMSKP'
tev_cat_res = [39, 74, 144] # H46, D81, C151


def get_distance(c1, c2):
	""" Returns the distance between two Rosetts XYZ coordinate vectors"""
	dist = sqrt((c2.x - c1.x) ** 2 + (c2.y - c1.y) ** 2 + (c2.z - c1.z) ** 2)
	return dist


def find_res_ca_coords(pose, resnum):
	""" For a given pose and residue number, returns the coordinates of CA """
	residue = pose.residue(resnum)
	CA = residue.atom('CA')
	return CA.xyz()


def find_cat_res(pose):
	""" 
	For a given pose, checks the coordinates of each CA for closest match to 
	the catalytic residues of HCV protease. Returns the corresponding residue 
	for each of the three members of the triad.
	"""
	HCV_coords = [find_res_ca_coords(hcv_pose, x) for x in cat_res]
	matching_residue_dists = ['none', 'none', 'none']
	matching_residue_numbers = ['none', 'none', 'none']
	chain_length = pose.total_residue()

	# Checking each residue against all members of the catalytic triad
	for resnum in range(1, chain_length + 1):
		res_coord = find_res_ca_coords(pose, resnum)
		for n, hcv_coord in enumerate(HCV_coords):
			distance = get_distance(hcv_coord, res_coord)
			if matching_residue_dists[n] > distance:
				matching_residue_dists[n] = distance
				matching_residue_numbers[n] = resnum

	# Listing matched residue numbers and residue types
	catalytic_matches = {}
	for res in matching_residue_numbers:
		catalytic_matches[res] = str(pose.residue(res).name1())

	return catalytic_matches


def get_secstruct(pose):
	""" Uses DSSP to get a secondary structure string for the given pose """
	sec_struct = Dssp(pose)
	sec_struct.insert_ss_into_pose(pose)
	ss = str(pose.secstruct())

	return ss


#def b_blocks(secstruct):
#	""" 
#	Reads through a given DSSP string and returns a list of b-sheet sections.
#	If a single residue is considered a loop between sheet residues, it is 
#	treated as part of the sheet. 
#	"""

def align_proteins(pose_1, pose_2):
	""" 
	Aligns two proteins usinf the Superimpose mover. The second protein will be 
	aligned to the first. Accommodates proteins of different sizes by aligning 
	only the number of residues in the smaller protein. This assumes that the 
	difference between the sizes is much smaller than the total number of 
	residues in each.
	"""
	# Determining the size of the smaller pose
	min_length = min(pose_1.total_residue(), pose_2.total_residue())

	# Creating mover from start to the smaller end, aligning with all BB atoms
	sim = SuperimposeMover(pose_1, min_length, 190, min_length, 190, False)

	sim.apply(pose_2)
	return


class protease_info():
	def __init__(self, qname='TEV', sname, qpose, spose):
		self.query_name = qname.upper()
		self.subject_name = sname.upper()
		self.query_pose = qpose
		self.subject_pose = spose

		# Dali info
		self.Z_score = None
		self.rmsd = None
		self.lali = None
		self.nres = None
		self.pID = None
		self.description = None

		# Alignment
		self.alignment = None
		self.aligned_residues = []

		# Catalytic triad
		self.nucleophile_res = None
		self.catalytic_his = None
		self.catalytic_acid = None

	def auto_calculate(self, dali_file='dali/0000_dali_pdb90_tev.txt', 
			align_file='dali/0000_seq_align.txt', structure_map=tev_map):
		""" Run all the calculation functions """
		self.get_dali_info(dali_file=dali_file)
		self.get_alignment(align_file=align_file)
		self.map_aligned_residues()
		self.map_structure_elements(structure_map=structure_map)

	def get_dali_info(self, dali_file='dali/0000_dali_pdb90_tev.txt'):
		"""
		Read in appropriate summary from Dali download about this protein, 
		including Z score (indicating structural similarity to the query 
		structure), RMSD to TEV protease (the original query), lali (the number 
		of structurally equivalent CA atoms), nres (the total number of 
		residues in the chain), pID (percentage of identical amino acids in 
		equivalent residues), and PDB description.
		Header line:
		Chain   Z    rmsd lali nres  %id Description
		"""
		# Read in Dali summary
		with open(dali_file, 'r') as r:
			match_summaries = r.readlines()

		# Find appropriate line in the summary by PDB name, stopping when found
		summary_line = None
		for ms in match_summaries:
			if self.subject_name in ms.upper:
				summary_line = ms.split()
				break

		# If no appropriate line is found, print error message and exit
		if summary_line == None:
			print("No matching protein identified in Dali summary")
			return

		# If line was found, read in its values
		dali_values = summary_line.split()
		self.Z_score = dali_values[1]
		self.rmsd = dali_values[2]
		self.lali = dali_values[3]
		self.nres = dali_values[4]
		self.pID = dali_values[5]
		self.description = ' '.join(dali_values[6:])

		return

	def get_alignment(self, align_file='dali/0000_seq_align.txt'):
		"""
		Read in sequence alignment file as a set of contiguous strings
		Hacky--may need tweaking to generalize.

		Alignment file has sets of five lines, with each set covering 60 
		residues in the alignment. The first line (0) is the secondary 
		structure of the query (TEV). The second (1) is the query sequence. The
		third (2) is the identity match (indicated as positive by |). The 
		fourth (3) is the subject sequence. The fifth (4) is the subject 
		secondary structure.
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
				if self.subject_name in sa.upper():
					begin_block = n

		# Extracting relevant text block
		alignment_block = seq_aligns[begin_block:end_block]
		# Cleaning block
		abclean = [i.strip() for i in alignment_block[2:] if i != '\n']
		# Chack that there are the right number of lines
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

		# Verifying all lines are equal length
		line_lengths = [len(i) for i in align_lines.values()]
		assert all([elem == line_lengths[0] for elem in line_lengths])

		self.alignment = list(align_lines.values())
		return

	def map_aligned_residues(self):
		"""
		Feed alignment data into a list of aligned_residues, each with 
		corresponding information about the position in both TEV and the 
		protease being analyzed.

		Start by finding the first five residues of the alignment in the pose 
		sequence for both TEV and the subject (adding 1 because Rosetta is 
		1-indexed). This is necessary because Dali only includes aligned 
		regions, whereas some subjects might have large N-terminal domains.
		"""
		# Get TEV residue number for beginning of alignment
		quer_match_seq = self.alignment[1].replace('-','').upper() 
		quer_pose_seq = self.query_pose.sequence()
		quer_pose_num = quer_pose_seq.find(quer_match_seq[:6]) + 1 

		# Get subject residue number for beginning of alignment
		subj_match_seq = self.alignment[3].replace('-','').upper() 
		subj_pose_seq = self.subject_pose.sequence()
		subj_pose_num = subj_pose_seq.find(subj_match_seq[:6]) + 1 

		# Loop through each residue in the alignment, adding aligned_residue
		# objects to the list for each
		for i in range(len(self.alignment[0])):
			# Getting PDB residue numbers
			quer_res_num = self.query_pose.residue(quer_pose_num).name1()
			subj_res_num = self.subject_pose.residue(quer_pose_num).name1()

			# Populating aligned_residue object
			aligned_residue = aligned_residue(
				quer_pose_num, quer_res_num, 
				subj_pose_num, subj_res_num, 
				*[row[i] for row in self.alignment])

			# Adding aligned_residue object to self.aligned_residues
			self.aligned_residues.append(aligned_residue)

			# Incrementing pose numbers
			quer_pose_num += 1
			subj_pose_num += 1

		return

	def map_structure_elements(self, structure_map=tev_map):


class aligned_residue():
	"""
	Data storage structure for a single residue. Includes information about 
	both the target residue in its own protein and the corresponding aligned 
	residue in the query structure. Information includes secondary structure, 
	whether residues are structurally matched (as opposed to unaligned), and 
	whether the residues are identical. Also stores residue numbers (both PDB
	and pose) for both residues.
	"""
	def __init__(self, qpnum, qrnum, spnum, srnum, qdss, qseq, qid, sseq, sdss):
		self.pose_number = spnum
		self.residue_number = srnum
		self.sec_struct = sdss
		self.res_type = sseq.upper()

		self.tev_pose_number = qpnum
		self.tev_residue_number = qrnum
		self.tev_sec_struct = qdss
		self.tev_res_type = qseq.upper()

		# Determine res identity, based on whether connection line was drawn
		if qid == '|':
			self.residues_equal = True
		else:
			self.residues_equal = False

		# Determine whether residues are structurally aligned, based on case
		if   all([i == i.upper() for i in [qseq, sseq]]):
			self.residues_align = True
		elif all([i == i.lower() for i in [qseq, sseq]]):
			self.residues_align = False
		else:
			print('Residue cases do not match')
			print(spnum, sseq, qpnum, qseq)
			assert False







