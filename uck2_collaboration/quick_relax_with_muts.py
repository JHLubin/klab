"""
cdkl; python uck2_collaboration/quick_relax_with_muts.py --ligand c35amc --mutant 39 Y --n 5
"""
import argparse
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	RestrictToRepackingRLT, RestrictToRepacking, \
	RestrictAbsentCanonicalAASRLT, OperateOnResidueSubset
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.select.residue_selector import \
	ResidueIndexSelector, NotResidueSelector, OrResidueSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.protocols.relax import FastRelax

parser = argparse.ArgumentParser()
parser.add_argument('--ligand', help='Which ligand?')
parser.add_argument('--mutant', nargs=2, 
	help='site residue (res 65 is 39 in pose)')
parser.add_argument('--n', type=int, help='How many decoys?')
args = parser.parse_args()

opts = '-use_input_sc -ex1 -ex2 -ignore_zero_occupancy false \
	-nblist_autoupdate -extra_res_fa \
	uck2_collaboration/params/{}.params'.format(args.ligand)
init(opts)

# Create score function
sf=create_score_function('ref2015_cst')
	# Upweight coordinate constraint because clashes are so bad
#sf.set_weight(ScoreType.coordinate_constraint, 2) 

# Make coordinate constraints mover
cg = CoordinateConstraintGenerator()
ac = AddConstraints()
ac.add_generator(cg)

# Make minimization mover, applied before FastRelax
minmov = MinMover()
minmov.min_type('lbfgs_armijo_nonmonotone')
	# Movemap will allow minimization of side chains only
mm = MoveMap()
mm.set_bb(False)
mm.set_chi(True)
minmov.movemap(mm)

# Make FastRelax mover
fr=FastRelax()
fr.set_scorefxn(sf)
	# Task factory for FastRelax
tf=TaskFactory()
if args.mutant:
	# Keep exclusion list of mutated residues
	not_mutated = NotResidueSelector()
	# If a mutations are input, force packing only to target residues
	res_selection = ResidueIndexSelector(str(args.mutant[0]))
	restriction = RestrictAbsentCanonicalAASRLT()
	restriction.aas_to_keep(args.mutant[1].upper())
	tf.push_back(OperateOnResidueSubset(restriction, res_selection))
	unmutated = NotResidueSelector(res_selection)
	tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), unmutated))
else:
	tf.push_back(RestrictToRepacking)

fr.set_task_factory(tf)

filnam = 'uck2_collaboration/start_models/{}.pdb'.format(args.ligand)
pose = pose_from_pdb(filnam)
ac.apply(pose)
minmov.apply(pose)

# Collect RMSD info in pose
rmm = RMSDMetric(pose)

jnam = '/scratch/jhl133/uck2_collaboration/check_muts/{}_{}{}'.format(args.ligand, *args.mutant)
jd = PyJobDistributor(jnam, args.n, sf)
while not jd.job_complete:
	pp = Pose(pose)
	fr.apply(pp)
	rmm.apply(pp)
	jd.output_decoy(pp)