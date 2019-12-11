from analyze_design_decoys import fix_file
import argparse
import design_protease as dp
from os import makedirs, remove
from os.path import basename, isdir, isfile, join
from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import \
	ChainSelector, OrResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--start_struct", required=True,
		default='fibrils_collaboration/other_models/crude_ext_cat.pdb', 
		help="Pick starting PDB")
	parser.add_argument("-od", "--out_dir", required=True,
		help="Name an output directory for decoys")
	parser.add_argument("-site", "--site", required=True, type=int,
		help="What is the asyn cleavage site in this frame?")
	parser.add_argument("-seq", "--sequence", required=True, type=str,
		help="What substrate sequence do you want to thread?")
	parser.add_argument("-n", "--number_decoys", type=int, default=10, 
		help="How many decoys should be made? (Default is 10.)")
	parser.add_argument('-x', '--extend', 
		help='Extend output name, ex: job from SLURM')
	args = parser.parse_args()
	assert len(args.sequence) == 12
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
	ors = OrResidueSelector()
	ors.add_residue_selector(ChainSelector('A')) # Constrain main backbone
	ors.add_residue_selector(ResidueIndexSelector('215-217')) # Preserving original peptide
	cg.set_residue_selector(ors)

	ac = AddConstraints()
	ac.add_generator(cg)
	ac.apply(pose)

	return 

args = parse_args()

if not isdir(args.out_dir):
	makedirs(args.out_dir)

opts = '-enzdes::cstfile fibrils_collaboration/htra1_protease.cst -run:preserve_header -mute core'
opts += ' -pep_refine -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false -score:weights ref2015_cst'
init(opts)

# Score function and starting PDB
sf = create_score_function('ref2015_cst')
pose = pose_from_pdb(args.start_struct)

# Setting FoldTree
ft=FoldTree()
ft.add_edge(1,211,-1)
ft.add_edge(1,217,1)
ft.add_edge(217,212,-1)
ft.add_edge(217,223,-1)
pose.fold_tree(ft)

# Changing peptide sequence
asyn_seq = args.sequence.upper()
pose = dp.make_residue_changes(pose, sf, asyn_seq, 212, [61, 91, 169], None)

# Creating FlexPepDock protocol using init options
fpdock = FlexPepDockingProtocol()

# Making name
decoy_name = join(args.out_dir, 'htra1_prot_asyn_ext')
decoy_name += '_' + str(args.site)
decoy_name += '_' + asyn_seq
if args.extend:
	decoy_name += '_' + args.extend

# Fixing constraints text block, since enzdes constraints are not dynamic
if asyn_seq[5] != 'T':
	fix_pdb = decoy_name + '.pdb'
	pose.dump_pdb(fix_pdb)
	fix_file(fix_pdb)
	pose = pose_from_pdb (fix_pdb)
	remove(fix_pdb)

# Applying constraints to the pose	
apply_constraints(pose)

jd = PyJobDistributor(decoy_name, args.number_decoys, sf)
while not jd.job_complete:
	pp = Pose(pose)
	fpdock.apply(pp)
	jd.output_decoy(pp)
