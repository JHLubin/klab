from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import \
	ChainSelector, OrResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
import design_protease as dp

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
	ors.add_residue_selector(ResidueIndexSelector('113-115')) # Preserving original peptide
	cg.set_residue_selector(ors)

	ac = AddConstraints()
	ac.add_generator(cg)
	ac.apply(pose)

	return 

opts = '-enzdes::cstfile fibrils_collaboration/htra1_pdz.cst -run:preserve_header -mute core'
opts += ' -pep_refine -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false -score:weights ref2015_cst'
init(opts)

# Score function and starting PDB
sf = create_score_function('ref2015_cst')
pose = pose_from_pdb('fibrils_collaboration/other_models/crude_ext_pdz.pdb')

# Setting FoldTree
ft=FoldTree()
ft.add_edge(1,105,-1)
ft.add_edge(1,115,1)
ft.add_edge(115,106,-1)
pose.fold_tree(ft)

# Changing peptide sequence
asyn_seq = "EGYQDYEPEA"
pose = dp.make_residue_changes(pose, sf, asyn_seq, 106, None, None)

# Applying constraints to the pose
apply_constraints(pose)

# Creating FlexPepDock protocol using init options
fpdock = FlexPepDockingProtocol()

decoy_name = 'fibrils_collaboration/relax_simulations/ext_asyn_pdz_20191209/htra1_pdz_asyn_ext_4'
jd = PyJobDistributor(decoy_name, 25, sf)
while not jd.job_complete:
	pp = Pose()
	pp.assign(pose)
	fpdock.apply(pp)
	jd.output_decoy(pp)
