import argparse
from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import InterGroupInterfaceByVectorSelector, ChainSelector, ResidueIndexSelector, OrResidueSelector, NotResidueSelector, AndResidueSelector
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
	PreventRepackingRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.teaching import SimpleThreadingMover
from pyrosetta.rosetta.protocols.relax import FastRelax

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

def shell_selection(selection_1, selection_2, base_range):
	shell_select = InterGroupInterfaceByVectorSelector()
	shell_select.group1_selector(selection_1) 
	shell_select.group2_selector(selection_2) 
	shell_select.nearby_atom_cut(base_range)
	shell_select.vector_dist_cut(base_range + 2)

	return AndResidueSelector(shell_select, NotResidueSelector(selection_2))

parser = argparse.ArgumentParser()
parser.add_argument('model', help='Select starting model.')
parser.add_argument('job', help='job from SLURM')
args = parser.parse_args()

#opts = '-use_input_sc -ex1 -ex2 -enzdes::cstfile htra1_pdz.cst -run:preserve_header'
opts = '-enzdes::cstfile htra1_pdz.cst -run:preserve_header'
init(opts)

pdz = ChainSelector('A')
peptide = ChainSelector('B')
designable = shell_selection(pdz, peptide, 6)
packable = shell_selection(pdz, designable, 6)
mobile = OrResidueSelector(designable, packable)
mobile.add_residue_selector(peptide)
static = NotResidueSelector(mobile)

sf = get_fa_scorefxn()

tf = TaskFactory()
tf.push_back(IncludeCurrent())
tf.push_back(ExtraRotamers(0, 1, 1))
tf.push_back(ExtraRotamers(0, 2, 1))

prevent = PreventRepackingRLT() # No repack, no design
repack = RestrictToRepackingRLT() # No design
tf.push_back(OperateOnResidueSubset(prevent, static))
tf.push_back(OperateOnResidueSubset(repack, packable))
tf.push_back(OperateOnResidueSubset(repack, peptide))

fr = FastRelax()
fr.set_scorefxn(sf)

fd = FastDesign()
fd.set_scorefxn(sf)
fd.set_task_factory(tf)

pose = pose_from_pdb(args.model)
pose.fold_tree(make_fold_tree())
pose = apply_constraints(pose)
pose = thread_to_htra1('QDYEPEA', pose)

jnam = 'pdz_designs/{}_designed_{}.pdb'.format(args.model, args.job)
fr.apply(pose)
fd.apply(pose)
pose.dump_pdb(jnam)
