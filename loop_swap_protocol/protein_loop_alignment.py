#!/usr/bin/python
from pyrosetta import *
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from math import sqrt

init()
hcv_pose = pose_from_pdb('a_to_s_ly104_WT.pdb')
cat_res = {72:'H', 96:'D', 154:'S'}

def get_distance(c1, c2):
	""" Returns the distance between two Rosetts XYZ coordinate vectors"""
	dist = sqrt((c2.x - c1.x) ** 2 + (c2.y - c1.y) ** 2 + (c2.z - c1.z) ** 2)
	return dist


def find_res_ca(pose, resnum):
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
	HCV_coords = [find_res_ca(hcv_pose, x) for x in cat_res]
	matching_residue_dists = ['none', 'none', 'none']
	matching_residue_numbers = ['none', 'none', 'none']
	chain_length = pose.total_residue()

	# Checking each residue against all members of the catalytic triad
	for resnum in range(1, chain_length + 1):
		res_coord = find_res_ca(pose, resnum)
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


def b_blocks(secstruct):
	""" 
	Reads through a given DSSP string and returns a list of b-sheet sections.
	If a single residue is considered a loop between sheet residues, it is 
	treated as part of the sheet. 
	"""