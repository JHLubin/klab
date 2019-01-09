"""
Using the amino acid sequence from a fasta for a fibrillary protein (such as 
tau, 5o3l), this script threads a sliding frame of residues from the fibril 
onto the substrate of the given PDB file (such as Htra1 protease, 3nzi) then 
runs a FastRelax. The backbone of the peptide is mobile, but had coordinate 
constraints applied. The aim is to determine the most favorable docking points 
along the chain, based only on sequence.
"""
import argparse
from os import makedirs
from os.path import basename, isdir, isfile, join
from pyrosetta import *
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pack.task.operation import \
	ExtraRotamers, IncludeCurrent, RestrictToRepacking
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.teaching import SimpleThreadingMover
from pyrosetta.toolbox import mutate_residue
from sys import exit

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--start_struct", required=True,
		help="Pick starting PDB")
	parser.add_argument("-o", "--out_dir", required=True,
		help="Name an output directory for decoys")
	parser.add_argument("-f", "--fasta", type=str, required=True,
		help="Pick a fasta sequence file to scan")
	parser.add_argument("-c", "--constraints", type=str, required=True,
		help="Pick a constraints file for the enzyme")
	parser.add_argument('-n', "--n_decoys", type=int, default=100, 
		help="How many decoys do you want? (Default: 100)")
	args = parser.parse_args()
	return args


def apply_constraints(pose):
	""" 
	Applies enzdes constraints form the input CST file to a pose
	Also applies coordinate constraints to the substrate peptide, assumed to 
	be chain B
	"""
	# Enzdes constraints
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)

	# Coordinate constraints
	cg = CoordinateConstraintGenerator()
	cs = ChainSelector('B')
	cg.set_residue_selector(cs)

	ac = AddConstraints()
	ac.add_generator(cg)
	ac.apply(pose)

	return 


def setup_fastrelax(sf, pose):
	"""
	Creates FastRelax mover with appropriate score function, movemap, and 
	packer rules. All sidechains are mobile, only substrate backbone is mobile.
	"""
	relax = FastRelax()
	relax.set_scorefxn(sf)

	# MoveMap
	mm = MoveMap()
	mm.set_chi(True)
	for i in range(1, pose.total_residue() + 1):
		if pose.residue(i).chain() == 2:
			mm.set_bb(i, True)
	relax.set_movemap(mm)

	# Packer tasks
	tf = standard_task_factory()
	tf.push_back(RestrictToRepacking())
	tf.push_back(IncludeCurrent())
	tf.push_back(ExtraRotamers(0, 1, 1))
	tf.push_back(ExtraRotamers(0, 2, 1))
	relax.set_task_factory(tf)

	return relax


def read_fasta(fasta):
	""" 
	Extracts the sequence from a fasta file, assumiong it's the second line 
	"""
	with open(fasta, 'r') as r:
		flines = r.readlines()

	return flines[1]


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
	opts = '-cst_fa_weight 1.0 -run:preserve_header -enzdes::cstfile {}'
	init(opts.format(args.constraints))

	# Score function and starting PDB
	sf = create_score_function('ref2015_cst')
	pose = pose_from_pdb(args.start_struct)

	# Applying constraints to the pose, if given
	apply_constraints(pose)

	# Making FastRelax mover
	fr = setup_fastrelax(sf, pose)

	# Reading in sequence
	scan_seq = read_fasta(args.fasta)

	# Going through all 6-residue frames in tau
	for frame in range(len(scan_seq))[:-4]:
		# Making name from the scissile bond and the frame sequence
		position_name = '{:03d}'.format(frame + 5)
		seq = scan_seq[frame:frame + 5]
		set_name = '_'.join([position_name, seq])
		print(set_name)

		# Make threaded model
		threaded_pose = thread_to_htra1(seq, pose)
		threaded_name = join(args.out_dir, 'threaded_htra1_' + set_name + '.pdb')
		if not isfile(threaded_name):
			threaded_pose.dump_pdb(threaded_name)

		# Creating relaxed decoys
		decoy_name = join(args.out_dir, 'htra1_scan_' + set_name)
		jd = PyJobDistributor(decoy_name, args.n_decoys, sf)
		while not jd.job_complete:
			pp = Pose()
			pp.assign(threaded_pose)
			fr.apply(pp)
			jd.output_decoy(pp)

if __name__ == '__main__':
	args = parse_args()
	main(args)