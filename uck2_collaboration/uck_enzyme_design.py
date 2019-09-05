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
from sys import exit

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--ligand', required=True, help='Which ligand?')
parser.add_argument('-od', '--outdir', type=str, required=True,help='Output directory?')
parser.add_argument('-ala', "--all_ala", action="store_true", help="Prepack to alanine?")
parser.add_argument('-x', '--extend', default='', help='Extend output name, ex: job from SLURM')
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
mutable_res = ResidueIndexSelector('39,42,47,53,55,57,62,63,86,91,111,112,140,144,150,158,163')
packable_res = ResidueIndexSelector('10,11,12,13,14,16,37,40,41,43,45,46,48,49,50,51,52,54,56,58,68,85,87,90,92,110,113,136,139,141,145,146,147,149,151,154,155,156,157,159,160,161,162,164,165,166,167,170,171')
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
for i in [1,2,3,4,5,6,7,8,9,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,38,44,59,60,61,62,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,88,89,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,142,143,148,152,153,168,169,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204]:
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

# Load and prepare pose
filnam = 'uck2_collaboration/start_models/{}.pdb'.format(args.ligand)
pose = pose_from_pdb(filnam)
# Pre-alanine packer
if args.all_ala:
	pt = tf_alla.create_task_and_apply_taskoperations(pose)
	aa_prm = PackRotamersMover(sf, pt)
	aa_prm.apply(pose)
minmov.apply(pose)

jnam = '{}/{}_designed_{}'.format(args.outdir, args.ligand, args.extend)
jd = PyJobDistributor(jnam, args.num_jobs, sf)
while not jd.job_complete:
	pp = Pose(pose)
	for i in range(4):
		dtask = tf_fr.create_task_and_apply_taskoperations(pp)
		fr_prm = PackRotamersMover(sf, dtask)
		fr_prm.apply(pp)
	jd.output_decoy(pp)

