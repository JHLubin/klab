from math import sqrt
from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import \
	AndResidueSelector, NeighborhoodResidueSelector, NotResidueSelector, \
	OrResidueSelector, ResidueIndexSelector, ResidueNameSelector 

def selector_to_list(pose, selector):
	selection_vector = selector.apply(pose)
	selection_list = []

	for i in range(len(selection_vector)):
		if selection_vector[i + 1] == 1:
			selection_list.append(i + 1)

	return selection_list


def isolate_main_chain(pose):
	""" 
	In case a substrate chain is present in the PDB, isolate only the main 
	chain. This function assumes that the largest chain in the pose is the 
	main chain; will not produce the desired result for multi-chain proteins. 

	Input: pose (full)
	Output: pose (only the largest chain)
	"""
	# list of separate poses for each chain in the original pose
	chains = pose.split_by_chain() 

	# default assumption is that the first chain is the main chain
	main_chain = 1 

	# check through all chains in the pose
	for i in range (1, pose.num_chains() + 1): 
		if len(pose.chain_sequence(i))>len(pose.chain_sequence(main_chain)): 
			main_chain = i

	return chains[main_chain]


def check_selector_nonzero_intersections(selectors_list, pose):
	"""
	Takes an input list of selectors and a pose, and checks the intersection
	of the selectors applied to the pose, returning a boolean indicating 
	whether the intersection set is empty, and the intersection set.
	"""
	# Make an AndResidueSelector
	ars = AndResidueSelector()

	# Add all selectors in the input list
	for sel in selectors_list:
		ars_s.add_residue_selector(sel)

	# Get intersection residues list
	intersection_list = selector_to_list(pose, ars)

	# Check whether list is empty
	is_intersection_empty = len(intersection_list) == 0

	return is_intersection_empty, intersection_list


def calculate_distance(point_1, point_2):
	""" 
	Calculates the distance between two points of any number of dimensions, 
	given as lists. Returns a float.
	"""
	# Verify same dimensions for both points, determine number
	assert len(point_1) == len(point_2)
	dimensions = len(point_1)

	# Calculate distance as sqrt(sum of (differences squared))
	difs_squared = [(point_1[i]-point_2[i]) ** 2 for i in range(dimensions)]
	dist = sqrt(sum(difs_squared))

	return dist


def determine_vector(point_1, point_2):
	""" 
	Determines a unit vector from two points, given as lists of coordinates 
	"""
	# Verify same dimensions for both points, determine number
	assert len(point_1) == len(point_2)
	dimensions = len(point_1)

	# Calculate difference between points
	difs = [(point_1[i]-point_2[i]) for i in range(dimensions)]

	# For vector length 1, divide components by distance between points
	norm = calculate_distance(point_1, point_2)

	return [difs[i] / norm for i in range(dimensions)]


def get_pose_atoms_distance(pose, res_1_num, res_1_atom, res_2_num, res_2_atom):
	"""
	Given a pose, two residue numbers, and a specified atom for each residue,
	return the distance between the specified atoms.

	Ex: get_atoms_distance(pose, 72, 'ND1', 154, 'CB')
	"""
	# Get coordinates of first atom
	r1 = pose.residue(res_1_num)
	a1 = r1.atom(r1.atom_index(res_1_atom))
	coords_1 = list(a1.xyz)

	# Get coordinates of second atom
	r2 = pose.residue(res_2_num)
	a2 = r1.atom(r2.atom_index(res_2_atom))
	coords_2 = list(a2.xyz)

	# Get distance between two points
	return calculate_distance(coords_1, coords_2)


def get_hist_ring_center(pose, hist_res):
	""" 
	Given a pose and a residue number for a histidine, determines the 
	coordinates for the geometric center of the 5-member ring
	"""
	res = pose.res(hist_res)
	
	# Make sure given resifue is a histidine
	assert res.name3() == 'HIS'

	# Determine coordinates of all members in the HIS ring
	atoms = ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
	coords = []
	for a in atoms:
		atom = res.atom(res.atom_index(a))
		coords.append(list(atom.xyz()))

	# Calculate midpoint as average of all coordinates
	center = [sum([p[i] for p in coords])/len(coords) for i in range(3)]
	return center


def get_cat_res(in_pose):
	"""
	For an input pose (presently assumed to be a serine protease), several 
	steps are performed to determine the catalytic triad residues.

	Firstly, the pose is stripped down to the largest chain, which should be 
	the main. This is to eliminate substrate chains that might falsly appear to 
	be part of the triad. 

	Secondly, selections of ser, his, and acid (asp and glu) residues are 
	collected. 

	Thirdly, a crude screen is conducted to weed out obviously false residues.
	The set of histidines is collected and scanned through individually. For 
	each histidine, the set of neighboring residues within 10 A is checked for 
	overlap with the sets of serines and acids, and if either set intersection 
	is empty, which means that there is no nearby serine and/or acid, then that 
	histidine cannot be part of the triad and is eliminated. Remaining 
	histidines are collected with the lists of their potential ser and acid 
	counterparts.

	Fourthly, the geometry of potential triads are examined. The atom distances 
	are checked between the ND's of the histidine and the respective CB of 
	serine and OD of the acid residue. The ND-CB distance should be 3.3 +/-
	0.5 A, and the ND-OD distance should be 2.6 +/- 0.5 A. 


	"""
	# 1. Isolate the protease main chain (no substrate)
	pose = isolate_main_chain(in_pose)

	# 2. Get selectors of serines, histidines, and acidic residues (aspartic and glutamic)
	target_res = ['SER', 'HIS', 'ASP,GLU']
	name_selectors = {}
	for tr in target_res:
		# Need to set name3 in selectors to get terminal residues, HIS_D, etc.
		rns = ResidueNameSelector()
		rns.set_residue_name3(tr)
		name_selectors[tr] = rns

	# 3. Scan through list of histidines to weed out any that don't have nearby 
		# serine and acid residues within 10 A
	potential_his = selector_to_list(pose, name_selectors['HIS'])
	potential_triads = {}
	for h in potential_his:
		# Make selector for all residues with CA within 10 A, excluding self
		nrs = NeighborhoodResidueSelector(ResidueIndexSelector(str(h)),10)
		nrs.set_include_focus_in_subset(False)

		# Select intersection of neighbors and serines
		no_match, potential_ser = \
			check_selector_nonzero_intersections([nrs, name_selectors['SER']], pose)
		if no_match:
			# If there are no nearby serines, this can't be catalytic hist
			continue

		# Select intersection of neighbors and acid residues
		no_match, potential_ac = \
			check_selector_nonzero_intersections([nrs, name_selectors['ASP,GLU']], pose)
		if no_match:
			# If there are no nearby acid residues, this can't be catalytic hist
			continue

		# Histidine has both requisite neighbors. Collect that residue and the 
			# lists of neighbors.
		potential_triads[h] = [potential_ser, potential_ac]

	# 4. Examine rotamer geometry of potential triads
	for h, [sers, acids] in potential_triads.items():
		for s in sers:
			for a in acids: 
				pass












