#!/usr/bin/python
"""
Using the tau amino acid sequence from a fibrillary protein (such as tau, 5o3l), 
this script threads a sliding frame of residues from the fibril onto the 
substrate of the given PDB file (such as Htra1 protease, 3nzi) then runs a 
FastRelax. The aim is to determine the most favorable docking points along the 
chain based only on sequence.
"""
import argparse
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

# Fibril protein sequences, with middle line of sequence comprising fibril core
tau_seq = 'AKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGK\
VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTF\
RENAKAKTDHGAEIVYKSPVV'

a_sys_seq = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGV\
LYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKK\
DQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'


def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--start_struct", required=True,
		help="Pick starting PDB")
	parser.add_argument("-o", "--out_dir", required=True,
		help="Name an output directory for decoys")
	parser.add_argument("-c", "--constraints", type=str, 
		help="Pick constraints file, if appropriate")
	parser.add_argument("-c", "--constraints", type=str, 
		help="Pick constraints file, if appropriate")


def apply_constraints(pose):
	""" Applies enzdes constraints form the input CST file to a pose """
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)
	return pose


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
	neighbors = [45, 46, 59, 60, 62, 124, 126, 145, 148, 149, 150, 157, 164, 
		165, 166, 167, 168, 170, 171, 183, 184, 185, 186, 187, 188, 189, 192, 
		193, 194, 195, 212, 213, 214, 215, 216] # Did 10A selection separately
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
	assert len(sequence) == 5

	# Constructing and applying mover
	tm = SimpleThreadingMover(sequence, 212)
	threaded_pose = Pose()
	threaded_pose.assign(pose)
	tm.apply(threaded_pose)

	return threaded_pose


def main(args):
	# Destination folder for PDB files
	if not isdir(args.out_dir):
		makedirs(args.out_dir)

	# Initialize Rosetta
	if args.constraints:
		opts = '-enzdes::cstfile {} -cst_fa_weight 1.0 -run:preserve_header'
		init(opts.format(args.constraints))
	else: 
		init()

	# Score function and starting PDB
	sf = create_score_function('ref2015_cst')
	pose = pose_from_pdb(args.start_struct)

	# Applying constraints to the pose, if given
	if args.constraints:
		pose = apply_constraints(pose)

	# Making FastRelax mover
	fr = setup_fastrelax(sf)

	# Going through all 6-residue frames in tau
	for frame in range(len(tau_seq))[:-4]:
		# Making name from position within tau sequence and the frame sequence
		position_name = '{:03d}'.format(frame + 275) # First res is V270
		# 275 selected so that name reflects which residue is downstream scissile
		seq = tau_seq[frame:frame + 5]
		set_name = '_'.join([position_name, seq])
		print set_name

		# Make threaded model
		threaded_pose = thread_to_htra1(seq, pose)
		threaded_name = 'threaded_htra1_tau_' + set_name + '.pdb'
		#########threaded_pose.dump_pdb(join(out_dir, threaded_name))
		###### Don't always overwrite threaded

		# Creating relaxed decoys
		decoy_name = join(out_dir, 'relaxed_htra1_tau_' + set_name)
		jd = PyJobDistributor(decoy_name, 20, sf)
		while not jd.job_complete:
			pp = Pose()
			pp.assign(threaded_pose)
			fr.apply(pp)
			jd.output_decoy(pp)

if __name__ == '__main__':
	args = parse_args()
	main(args)