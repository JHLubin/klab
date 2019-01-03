#!/usr/bin/python
"""
Using the tau amino acid sequence from 5o3l, this script threads a sliding 
frame of 6 tau residues into the substrate of Htra1 protease 3nzi then runs a 
FastRelax. The aim is to determine the most favorable docking points along the
tau chain based only on sequence.
"""
from os import makedirs
from os.path import basename, isdir, isfile, join
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
	Make a fold tree that connects the res 385 to the 
	substrate terminal residue. More efficient sampling.
	Presently hard-coded for Htra1 PDZ
		I385 is residue 10.
		Last residue of PDZ chain A is 105
		Terminal residue of substrate chain B is 112
		Substrate chain B is 7 residues 
	"""
	ft = FoldTree()
	ft.add_edge(10, 1, -1)
	ft.add_edge(10, 105, -1)
	ft.add_edge(10, 112, 1)
	ft.add_edge(112 ,106, -1)
	assert ft.check_fold_tree()

	return ft


def setup_fastrelax(sf):
	"""
	Creates FastRelax mover with appropriate score function, movemap, and 
	packer rules. List of neighbor residues was generated using a 8A 
	PyMOL selection expansion around the peptide chain.
	"""
	relax = FastRelax()
	relax.set_scorefxn(sf)

	# MoveMap
	mm = MoveMap()
	mm.set_bb_true_range(106, 112)
	neighbors = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 
		22, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 53, 59, 62, 69, 70, 
		71, 72, 73, 74, 75, 76, 77, 78, 79, 83, 85, 98, 100, 101, 103, 104, 
		105, 106, 107, 108, 109, 110, 111, 112] # Did 8A selection separately
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
	assert len(sequence) == 7

	# Constructing and applying mover
	tm = SimpleThreadingMover(sequence, 106)
	threaded_pose = Pose()
	threaded_pose.assign(pose)
	tm.apply(threaded_pose)

	return threaded_pose


def main():
	# Destination folder for PDB files
	out_dir = 'pdz_tau_slide'
	if not isdir(out_dir):
		makedirs(out_dir)

	# Initialize Rosetta
	opts = '-enzdes::cstfile htra1_pdz.cst \
			-cst_fa_weight 1.0 -run:preserve_header'
	init(opts)

	# Score function
	sf = create_score_function('ref2015_cst')

	# Starting poses
	poses = []
	pdbs = ['pdz_relaxed_1.pdb', 'pdz_relaxed_2.pdb', 
			'pdz_relaxed_3.pdb', 'pdz_relaxed_4.pdb']
	for pdb in pdbs:
		p = pose_from_pdb(pdb)

		# Applying fold tree and constraints to the pose
		p.fold_tree(make_fold_tree())
		p = apply_constraints(p)

		poses.append(p) 

	# Making FastRelax mover
	fr = setup_fastrelax(sf)

	# Going through all 7-residue frames in tau
	for frame in range(len(tau_seq))[:-6]:
		# Making name from position within tau sequence and the frame sequence
		position_name = '{:03d}'.format(frame + 276) # First res is V270
		# P276 is the first frame's C-terminal residue
		seq = tau_seq[frame:frame + 7]
		set_name = '_'.join([position_name, seq])
		print set_name

		for n, pose in enumerate(poses):
			# Creating relaxed decoys
			r_name = 'relaxed_pdz_' + str(n) + '_tau_' + set_name
			jd = PyJobDistributor(join(out_dir, r_name), 50, sf)
			while not jd.job_complete:
				# Make threaded, relaxed models
				threaded_pose = thread_to_htra1(seq, pose)
				fr.apply(threaded_pose)
				jd.output_decoy(threaded_pose)

if __name__ == '__main__':
	main()