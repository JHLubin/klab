import argparse
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	ExtraRotamers, IncludeCurrent, OperateOnResidueSubset, \
	PreventRepackingRLT, RestrictAbsentCanonicalAASRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.core.select.residue_selector import \
	OrResidueSelector, NeighborhoodResidueSelector, NotResidueSelector, \
	ResidueIndexSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.protocols.minimization_packing import MinMover, \
	PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from sys import exit

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--ligand', required=True, help='Which ligand?')
parser.add_argument('-od', '--outdir', type=str, required=True,help='Output directory?')
parser.add_argument('-ala', "--all_ala", action="store_true", help="Prepack to alanine?")
parser.add_argument('-x', '--extend', help='Estend output name, ex: job from SLURM')
parser.add_argument('-n', '--num_jobs', type=int, default=1, help='How many decoys?')
args = parser.parse_args()

# Initialize
opts = '-use_input_sc -ex1 -ex2 -ignore_zero_occupancy false \
	-nblist_autoupdate -extra_res_fa \
	uck2_collaboration/params/{}.params'.format(args.ligand) 
init(opts)

# Selections based on PyMOL manual picks
ligand = ResidueIndexSelector(205)
catalytic_res = ResidueIndexSelector('15,33,143,148')
#mutable_res = ResidueIndexSelector('35,39,62,63,88')
mutable_res = ResidueIndexSelector('11,12,37,39,42,46,47,50,53,54,55,56,57,58,86,91,110,111,112,140,143,144,147,149,150,155,158,159,162,163,166')
#packable_res = ResidueIndexSelector('34,36,37,38,40,41,55,56,57,58,59,60,61,64,65,68,86,87,89,90,91,93,108,109,110,111,112,113,148')
packable_res = ResidueIndexSelector('10,13,14,16,35,40,41,43,45,48,49,51,52,54,32,85,87,90,92,113,136,139,141,145,146,151,154,156,157,160,161,164,165,167,170,171')
mobile_res = OrResidueSelector(mutable_res, packable_res)
mobile_res.add_residue_selector(ligand)
static_res = NotResidueSelector(mobile_res)

# Make minimization mover, applied before FastRelax
sf = get_fa_scorefxn()
minmov = MinMover()
minmov.score_function(sf)
minmov.min_type('lbfgs_armijo_nonmonotone')
	# Movemap will allow minimization of side chains only
mm = MoveMap()
mm.set_bb(False)
mm.set_chi(True)
for i in [15,33,143,148]:
	mm.set_chi(i, False)
minmov.movemap(mm)

# Task factory for FastRelax
prevent = PreventRepackingRLT() # No repack, no design
repack = RestrictToRepackingRLT() # No design
tf_fr = TaskFactory()
tf_fr.push_back(IncludeCurrent())
tf_fr.push_back(ExtraRotamers(0, 1, 1))
tf_fr.push_back(ExtraRotamers(0, 2, 1))
tf_fr.push_back(OperateOnResidueSubset(repack, packable_res))
tf_fr.push_back(OperateOnResidueSubset(repack, ligand))
tf_fr.push_back(OperateOnResidueSubset(prevent, static_res))

if args.all_ala:
	# Task factory for pre-alanine
	tf_alla = TaskFactory()
	restriction = RestrictAbsentCanonicalAASRLT()
	restriction.aas_to_keep('A')
	tf_alla.push_back(OperateOnResidueSubset(restriction, mutable_res))
	tf_alla.push_back(OperateOnResidueSubset(repack, packable_res))
	tf_alla.push_back(OperateOnResidueSubset(repack, ligand))
	tf_alla.push_back(OperateOnResidueSubset(prevent, static_res))

# Set up FastRelax
fr = FastRelax()
fr.set_scorefxn(sf)
fr.set_movemap(mm)
fr.set_task_factory(tf_fr)

# Load and prepare pose
filnam = 'uck2_collaboration/start_models/{}.pdb'.format(args.ligand)
pose = pose_from_pdb(filnam)
# Pre-alanine packer
if args.all_ala:
	pt = tf_alla.create_task_and_apply_taskoperations(pose)
	prm = PackRotamersMover(sf, pt)
	prm.apply(pose)
minmov.apply(pose)

jnam = '{}/{}_designed_0{}'.format(args.outdir, args.ligand, args.extend)
jd = PyJobDistributor(jnam, args.num_jobs, sf)
while not jd.job_complete:
	pp = Pose(pose)
	fr.apply(pp)
	jd.output_decoy(pp)

