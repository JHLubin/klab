from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import \
	ChainSelector, OrResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol

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
	ors.add_residue_selector(ResidueIndexSelector('215-217')) # Preserving b-sheet templating
	cg.set_residue_selector(ors)

	ac = AddConstraints()
	ac.add_generator(cg)
	ac.apply(pose)

	return 

opts = '-enzdes::cstfile htra1_protease.cst -run:preserve_header'
opts += ' -pep_refine -ex1 -ex2 -use_input_sc -score:weights ref2015_cst'
init(opts)

# Score function and starting PDB
sf = create_score_function('ref2015_cst')
pose = pose_from_pdb('crude_ext_cat.pdb')

# Applying constraints to the pose
apply_constraints(pose)

# Creating FlexPepDock protocol using init options
fpdock = FlexPepDockingProtocol()

decoy_name = 'ext_cat_flexpep/htra1_protease_ext_subst_4'
jd = PyJobDistributor(decoy_name, 40, sf)
while not jd.job_complete:
	pp = Pose()
	pp.assign(pose)
	fpdock.apply(pp)
	jd.output_decoy(pp)
