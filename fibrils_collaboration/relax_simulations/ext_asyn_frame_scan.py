import argparse
from os.path import join
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task.operation import \
	ExtraRotamers, IncludeCurrent, RestrictToRepacking
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.teaching import SimpleThreadingMover

start_pdbs = ['htra1_protease_ext_subst_08_28.pdb', 
			'htra1_protease_ext_subst_03_47.pdb', 
			'htra1_protease_ext_subst_02_34.pdb', 
			'htra1_protease_ext_subst_03_16.pdb', 
			'htra1_protease_ext_subst_08_99.pdb', 
			'htra1_protease_ext_subst_08_66.pdb', 
			'htra1_protease_ext_subst_09_58.pdb']

pdb_dir = 'ext_cat_flexpep_bests'

# 12-res frames with cut sites 86-92
frames = ['TVEGAGSIAAAT',
		'VEGAGSIAAATG',
		'EGAGSIAAATGF',
		'GAGSIAAATGFV',
		'AGSIAAATGFVK',
		'GSIAAATGFVKK',
		'SIAAATGFVKKD']

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-p", "--pdb_number", required=True, type=int,
		help="Pick starting PDB from set")
	parser.add_argument("-f", "--frame_number", required=True, type=int,
		help="Pick frame to test")
	parser.add_argument('-n', "--name", type=str,
		help="Name variant")
	args = parser.parse_args()
	return args

def apply_constraints(pose):
	""" 
	Applies enzdes constraints from the input CST file to a pose
	Also applies coordinate constraints to the original substrate peptide
	which is assumed correct from the PDB
	"""
	# Enzdes constraints
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)

	# Coordinate constraints
	cg = CoordinateConstraintGenerator()
	cs = ResidueIndexSelector('213-216')
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


def thread_to_htra1(sequence, pose):
	"""
	Uses SimpleThreadingMover to swap out the native substrate and put in a
	new test sequence docked with Htra1 protease. The substrate peptide begins 
	at residue 212 of the pose, based on 3nzi.
	"""
	# Constructing and applying mover
	tm = SimpleThreadingMover(sequence, 212)
	threaded_pose = Pose()
	threaded_pose.assign(pose)
	tm.apply(threaded_pose)

	return threaded_pose



args = parse_args()

# Initialize Rosetta
opts = '-cst_fa_weight 1.0 -run:preserve_header -enzdes::cstfile htra1_protease.cst -score:weights ref2015_cst'
opts += ' -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false -cst_fa_weight 1.0 -pep_refine'
init(opts)

# Score function and starting PDB
sf = create_score_function('ref2015_cst')
pdb = join(pdb_dir,start_pdbs[args.pdb_number])
pose = pose_from_pdb(pdb)
apply_constraints(pose)

# Making FastRelax mover
fr = setup_fastrelax(sf, pose)

# Threading
seq = frames[args.frame_number]
threaded_pose = thread_to_htra1(seq, pose)

# Naming decoy
out_name = start_pdbs[args.pdb_number].replace('.pdb', '')
out_name += '_' + seq
out_name = join('ext_pep_asym_cut_check', out_name)
if args.name:
	out_name += args.name

# Relax threaded structure, run FlexPepDock
jd = PyJobDistributor(out_name, 100, sf)
while not jd.job_complete:
	pp = Pose()
	pp.assign(threaded_pose)
	fr.apply(pp)
	FlexPepDockingProtocol().apply(pp)
	jd.output_decoy(pp)