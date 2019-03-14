#!/usr/bin/python
"""
Testing the effects of pulling out parts of the tau fibril to get a sense of 
the energy cost, which might help get a sense of where the PDZ of Htra1 binds.
The approach is to pick a residue in the middle strand of the PDB model of the 
tau fibril, and apply constraints to push it 10A away from the corresponding 
residues in the chains above and below it and FastRelax-ing. This should pull 
a small stretch of the fibril out. Tugging will be repeated for every fourth 
residue, and 10 decoys will be produced.
"""
import argparse
from math import acos, radians
from os import makedirs
from os.path import isdir, join
from pyrosetta import *
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	OperateOnResidueSubset, PreventRepackingRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.constraints import \
	AtomPairConstraint, DihedralConstraint
from pyrosetta.rosetta.core.scoring.func import \
	CircularHarmonicFunc, HarmonicFunc
from pyrosetta.rosetta.core.select.residue_selector import \
	ChainSelector, NeighborhoodResidueSelector, \
	NotResidueSelector, OrResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.protocols.relax import FastRelax
from sys import exit

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--start_struct",
		default='tau_relaxed.pdb', help="Pick starting PDB")
	parser.add_argument("-o", "--outdir", 
		default='fibril_tugs', help="Name an output directory for decoys")
	parser.add_argument("-m", "--mobile_range", type=int, default=False, 
		help="Allow movement only X residues from tugged residue (ex: 5)")
	parser.add_argument("-d", "--tug_distance", type=float, default=6.5, 
		help="How far from the fibril the target residue should be tugged \
		(default: 6.5)")
	parser.add_argument("-n", "--num_decoys", type=int, default=50, 
		help="Desired number of decoys (default: 50)")
	group = parser.add_mutually_exclusive_group()
	group.add_argument("-t", "--test_single", type=int, default=False, 
		help="Pick a single residue to test tug on (ex: 347)")
	group.add_argument("-f", "--frame", type=int, default=1, 
		help="Tug one residue in per X residues (default: 1)")
	
	args = parser.parse_args()
	return args


def selector_to_list(pose, selector):
	""" Converts a selector output vector to a list of selected residues """
	selection_vector = selector.apply(pose)
	selection_list = []
	for i in range(len(selection_vector)): 
		if selection_vector[i+1]==1:
			selection_list.append(i+1)

	return selection_list 


def get_mobile_range(chain_no, central_res=None, range_from_tug=None):
	""" 
	Returns a selection of residues in given chain. The list can either include 
	the full chain, or a subset within a given range from a given central 
	residue. Still includes stom hard-coding for tau in the case of the subset.
	"""
	# Full chain
	if not all([central_res, range_from_tug]):
		chain_letter = chr(64 + chain_no) # 'A' is chr 65, 'B' is 66, etc.
		chain_selection = ChainSelector(chain_letter)
		return chain_selection

	# Subset of chain
	elif all([central_res, range_from_tug]):
		# Ignoring upstream residues of central res is near N-term
		if central_res >= 306 + range_from_tug:
			run_start = range_from_tug
		else:
			run_start = central_res - 306

		# Ignoring downstream residues if central res is near C-term
		if central_res <= 378 - range_from_tug:
			run_end = range_from_tug
		else:
			run_end = 378 - central_res

		run_str = '-'.join([str(site - run_start), str(site + run_end)])
		chain_run_selection = ResidueIndexSelector(run_str)	
		return chain_run_selection

	# Shouldn't happen, but just in case...
	else:
		print("Bug detected: received one input, not zero or two.")
		exit()


def get_bordering_atom_id(pose, site_chain, site_res, stack_position, atom):
	"""
	Determines the residue above or below the target residue in the tau fibril 
	and returns the specified atom of that neighbor. Selection assumes the  
	chain numbering matches 5o3l, wherein alternating chains are in each stack. 
	(Target stack is A-C-E-G-I.)
	"""
	# Verifying correct atom name
	if atom not in ['N', 'CA', 'C']:
		print("Incorrect atom name")
		exit()

	neighbor_chain_number = site_chain + 2 * stack_position
	neighbor_chain = chr(64 + neighbor_chain_number) # 'A' is chr(65)
	neighbor_res = pose.pdb_info().pdb2pose(neighbor_chain, site_res)
	atom_id = AtomID(pose.residue(neighbor_res).atom_index(atom), neighbor_res)

	return atom_id


def make_constraints(pose, chain_no, tug_site, tug_len):
	"""
	Creates harmonic constraints for 5o3l (hard-coded) to pull a target site 
	residue away from its neighbors in the chains above and below by a given 
	length in Angstroms. Uses get_bordering_atom_id to determine neighbors.
	"""
	# Catch useless tugs
	if tug_len < 5:
		print("Default distance is roughly 5A, so pick something bigger for d")
		exit()

	# Identifying target residue's CA
	res_info = pose.pdb_info().pose2pdb(tug_site).split()
	res_no_in_chain = int(res_info[0])
	ca_target = AtomID(pose.residue(tug_site).atom_index('CA'), tug_site)

	# Determining residues immediately above and below the target
	ca_ab = get_bordering_atom_id(pose, chain_no, res_no_in_chain, -1, 'CA')
	ca_bl =	get_bordering_atom_id(pose, chain_no, res_no_in_chain, 1, 'CA')

	# Creating repulsion constraints to push away from neighbors
	tug_sd = float(tug_len - 4.5) / 10
	hf = HarmonicFunc(tug_len, tug_sd) # Normally distance is <4.8A
	cst_above = AtomPairConstraint(ca_target, ca_ab, hf)
	cst_below = AtomPairConstraint(ca_target, ca_bl, hf)

	# Determining CA above and below target one res back for dihedrals
	ca_ab_p = get_bordering_atom_id(pose, chain_no, res_no_in_chain - 1, -1, 'CA')
	ca_bl_p = get_bordering_atom_id(pose, chain_no, res_no_in_chain - 1, 1, 'CA')

	# Making constraints--30 degree tolerance from ideal angle
	chf_above = CircularHarmonicFunc(radians(90), .5, .5)
	chf_below = CircularHarmonicFunc(radians(-90), .5, .5)

	dih_above = DihedralConstraint(ca_ab_p, ca_ab, ca_bl, ca_target, chf_above)
	dih_below = DihedralConstraint(ca_bl_p, ca_bl, ca_ab, ca_target, chf_below)

	return [cst_above, cst_below, dih_above, dih_below]


def make_move_map(pose, bb_selection, sc_selection):
	""" 
	Makes a movemap from a pose and two residue selections. All backbones are 
	fixed except those in the first selection, and only residues in the second 
	selection repackable, with other side chains fixed.
	"""
	mm = MoveMap()

	# Setting backbones
	for i in selector_to_list(pose, bb_selection):
		mm.set_bb(i, True)

	# Setting repacking
	for i in selector_to_list(pose, sc_selection):
		mm.set_chi(i, True)

	return mm


def make_task_factory(repackable_selection):
	""" 
	Creates a task factory with residues in a selection repackable (not 
	designable), and all other residues fixed.
	"""
	# Repack options
	prevent = PreventRepackingRLT() # No repack, no design
	repack = RestrictToRepackingRLT() # No design

	# Creating compliment selection
	fixed_selection = NotResidueSelector(repackable_selection)

	# Making task factory
	tf = TaskFactory()
	tf.push_back(OperateOnResidueSubset(prevent, fixed_selection))
	tf.push_back(OperateOnResidueSubset(repack, repackable_selection))

	return tf


def main(args):
	init()
	sf = create_score_function('ref2015_cst')
	sf.set_weight(ScoreType(1).atom_pair_constraint, 2) # Increasing repulsion weight
	sf.set_weight(ScoreType(1).dihedral_constraint, 0.5) # Reducing dihedral weight
	def_score = get_fa_scorefxn()

	if not isdir(args.outdir):
		makedirs(args.outdir)

	# Reading in PDB, determining middle chain
	pose = pose_from_pdb(args.start_struct)
	chain_count = pose.num_chains()
	chain_no = chain_count / 2

	# Determining which sites to tug
	if args.test_single:
		tug_sites = [args.test_single]
	else:
		chain_residues = selector_to_list(pose, get_mobile_range(chain_no))
		tug_sites = chain_residues[1::args.frame] # First res will error
		# Error because dihedrals are calculated using the upstream residue

	# Making decoy set for each sliding frame on the tau middle monomer
	for site in tug_sites:
		in_res = int(pose.pdb_info().pose2pdb(site).split()[0])
		site_name = join(args.outdir, 'tau_fibril_distort_' + str(in_res))

		# Determining backbone-mobile section of pose
		if args.mobile_range:
			mobile_selection = get_mobile_range(chain_no, \
				central_res=in_res, range_from_tug=args.mobile_range)
		else:
			mobile_selection = get_mobile_range(chain_no)

		# Selecting repackable shell within 18A of the target (includes target)
		ris = ResidueIndexSelector(str(site))
		nrs = NeighborhoodResidueSelector(ris, 18)

		# Combining selections
		combined_selection = OrResidueSelector()
		combined_selection.add_residue_selector(mobile_selection)
		combined_selection.add_residue_selector(nrs)

		# Making move map and task factory for FastRelax
		mm = make_move_map(pose, mobile_selection, combined_selection)
		tf = make_task_factory(combined_selection)

		# Making FastRelax
		fr = FastRelax()
		fr.set_scorefxn(sf)
		fr.set_movemap(mm)
		fr.set_task_factory(tf)

		# Making constraints
		csts = make_constraints(pose, chain_no, site, args.tug_distance)

		# Making ten decoys
		jd = PyJobDistributor(site_name, args.num_decoys, sf)
		while not jd.job_complete:
			p = Pose()
			p.assign(pose)

			# Add constraints
			for c in csts:
				p.add_constraint(c)

			fr.apply(p)

			print("\n" * 5, jd.current_name)
			print(sf.show(p), "\n" * 5)

			jd.output_decoy(p)

if __name__ == '__main__':
	args = parse_args()
	main(args)