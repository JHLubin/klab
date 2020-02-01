from math import sqrt
from pyrosetta import *
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.scoring import rmsd_atoms
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.protocols.simple_moves import SuperimposeMover

def get_secstruct(pose):
	""" Uses DSSP to get a secondary structure string for the given pose """
	sec_struct = Dssp(pose)
	sec_struct.insert_ss_into_pose(pose)
	ss = str(pose.secstruct())

	return ss


def get_section_RMSD(pose_1, selector_1, pose_2, selector_2):
	"""
	Calculates CA RMSD of pose regions. If selectors are given as none, will 
	calculate the whole pose.
	"""
	rmsd = RMSDMetric()
	rmsd.set_rmsd_type(rmsd_atoms.rmsd_protein_bb_ca)
	rmsd.set_comparison_pose(pose_1)
	if selector_1:
		rmsd.set_residue_selector_reference(selector_1)
	if selector_2:
		rmsd.set_residue_selector(selector_2)
	
	return rmsd.calculate(pose_2)


def make_subpose(pose, start=1, end=1):
	""" Return pose that is a section of another pose """
	# Set end if not given
	if end == 1:
		end = pose.total_residue()

	# Make new subpose
	subpose = Pose(pose, start, end)

	return subpose


def align_proteins(pose_1, pose_2):
	""" 
	Aligns two proteins using the Superimpose mover. The second protein will be 
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


#def align_proteins(pose_1, pose_2):
	""" 
	Aligns two proteins usinf the Superimpose mover. The second protein will be 
	aligned to the first. Accommodates proteins of different sizes by aligning 
	only the number of residues in the smaller protein. This assumes that the 
	difference between the sizes is much smaller than the total number of 
	residues in each.
	"""
	# Determining the size of the smaller pose
	#min_length = min(pose_1.total_residue(), pose_2.total_residue())

	# Creating mover from start to the smaller end, aligning with all BB atoms
	#sim = SuperimposeMover(pose_1, min_length, 190, min_length, 190, False)

	#sim.apply(pose_2)
	#return


#def get_rmsd(pose_1, pose_2, residues_1, residues_2):
	#assert len(residues_1) == len(residues_2)
	#n = len(residues_1)

	#difs = 0
	#for i in range(n):
	#	r1_coords = find_res_ca_coords(pose_1, residues_1[i])
	#	r2_coords = find_res_ca_coords(pose_2, residues_2[i])
	#	difs += (r2_coords.x - r1_coords.x) ** 2
	#	difs += (r2_coords.y - r1_coords.y) ** 2
	#	difs += (r2_coords.z - r1_coords.z) ** 2

	#return sqrt(difs / n)