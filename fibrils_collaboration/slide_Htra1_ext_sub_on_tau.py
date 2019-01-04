#!/usr/bin/python
"""
Using the tau amino acid sequence from 5o3l, this script threads a sliding 
frame of 6 tau residues into the substrate of Htra1 protease 3nzi then runs a 
FastRelax. The aim is to determine the most favorable docking points along the
tau chain based only on sequence.
"""
from os import makedirs
from os.path import basename, isdir, join
from pyrosetta import *
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pack.task.operation import \
	ExtraRotamers, IncludeCurrent, RestrictToRepacking
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.teaching import SimpleThreadingMover
from pyrosetta.toolbox import mutate_residue
from sys import exit

tau_seq = 'AKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGK\
VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTF\
RENAKAKTDHGAEIVYKSPVV'

def apply_constraints(pose):
	""" Applies enzdes constraints form the input CST file to a pose """
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)
	return pose


def make_fold_tree():
	"""
	Make a fold tree that connects the first catalytic serine to the 
	substrate scissile residue. More efficient sampling.
	Presently hard-coded for Htra1 protease
		S328 is residue 169.
		Last residue of protease chain A is 212
		Scissile residue of substrate chain B is 216
		Substrate chain B is residues 
	"""
	ft = FoldTree()
	ft.add_edge(169, 1, -1)
	ft.add_edge(169, 211, -1)
	ft.add_edge(169, 216, 1)
	ft.add_edge(216 ,212, -1)
	ft.add_edge(216 ,220, -1)
	assert ft.check_fold_tree()

	return ft


def setup_fastrelax(sf):
	"""
	Creates FastRelax mover with appropriate score function, movemap, and 
	packer rules. List of neighbor residues was generated using a 10A 
	neighborhood residue selector around the peptide chain.
	"""
	relax = FastRelax()
	relax.set_scorefxn(sf)

	# MoveMap
	mm = MoveMap()
	mm.set_bb_true_range(212,216)
	neighbors = [28, 29, 41, 42, 43, 44, 45, 46, 59, 60, 61, 62, 91, 125, 126, 
		127, 128, 129, 145, 147, 148, 149, 150, 151, 157, 164, 165, 166, 167, 
		168, 169, 170, 171, 183, 184, 185, 186, 187, 188, 189, 192, 193, 194, 212, 
		213, 214, 215, 216, 217, 218, 219, 220] # Did 10A selection separately
	for n in neighbors:
		mm.set_chi(n, True)
	relax.set_movemap(mm)

	# Packer tasks
	tf = standard_task_factory()
	tf.push_back(RestrictToRepacking())
	tf.push_back(IncludeCurrent())
	tf.push_back(ExtraRotamers(0, 1, 1))
	tf.push_back(ExtraRotamers(0, 2, 1))
	relax.set_task_factory(tf)

	return relax


def thread_to_htra1(sequence, pose):
	"""
	Uses SimpleThreadingMover to swap out the native substrate and put in a
	new test sequence docked with Htra1 protease. The substrate peptide begins 
	at residue 212 of the pose, based on 3nzi.
	"""
	assert len(sequence) == 9

	# Constructing and applying mover
	tm = SimpleThreadingMover(sequence, 212)
	threaded_pose = Pose()
	threaded_pose.assign(pose)
	tm.apply(threaded_pose)

	return threaded_pose


def main():
	# Destination folder for PDB files
	out_dir = 'ext_cat_tau_slide'
	if not isdir(out_dir):
		makedirs(out_dir)

	# Initialize Rosetta
	opts = '-enzdes::cstfile htra1_cat_general.cst \
			-cst_fa_weight 1.0 -run:preserve_header'
	init(opts)

	# Score function and starting PDB
	sf = create_score_function('ref2015_cst')
	pose = pose_from_pdb('ext_cat_relax.pdb')

	# Applying fold tree and constraints to the pose, deactivating by mutation
	pose.fold_tree(make_fold_tree())
	pose = apply_constraints(pose)
	mutate_residue(pose, 169, 'A') # Catalytic S328 (169 in pose) mutated to A 

	# Making FastRelax mover
	fr = setup_fastrelax(sf)

	# Going through all 6-residue frames in tau
	for frame in range(len(tau_seq))[:-8]:
		# Making name from position within tau sequence and the frame sequence
		position_name = '{:03d}'.format(frame + 279) # First res is V270
		# 275 selected so that name reflects which residue is downstream scissile
		seq = tau_seq[frame:frame + 9]
		set_name = '_'.join([position_name, seq])
		print(set_name)

		# Make threaded model
		threaded_pose = thread_to_htra1(seq, pose)
		threaded_name = 'threaded_htra1_tau_' + set_name + '.pdb'
		threaded_pose.dump_pdb(join(out_dir, threaded_name))

		# Creating relaxed decoys
		decoy_name = join(out_dir, 'relaxed_htra1_tau_' + set_name)
		jd = PyJobDistributor(decoy_name, 50, sf)
		while not jd.job_complete:
			pp = Pose()
			pp.assign(threaded_pose)
			fr.apply(pp)
			jd.output_decoy(pp)

if __name__ == '__main__':
	main()