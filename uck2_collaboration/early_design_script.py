import argparse
from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, OrResidueSelector, NotResidueSelector
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
	PreventRepackingRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign

parser = argparse.ArgumentParser()
parser.add_argument('ligand', help='Which ligand?')
parser.add_argument('job', help='job from SLURM')
args = parser.parse_args()

opts = '-use_input_sc -ex1 -ex2 -ignore_zero_occupancy false -nblist_autoupdate  -extra_res_fa '
opts += 'params/{}.params'.format(args.ligand) 
init(opts)

mutable_res = ResidueIndexSelector('35,39,62,63,88')
packable_res = ResidueIndexSelector('34,36,37,38,40,41,55,56,57,58,59,60,61,64,65,68,86,87,89,90,91,93,108,109,110,111,112,113,148')
ligand = ResidueIndexSelector(205)
mobile_res = OrResidueSelector(mutable_res, packable_res)
mobile_res.add_residue_selector(ligand)
static_res = NotResidueSelector(mobile_res)

sf = get_fa_scorefxn()

tf = TaskFactory()
tf.push_back(IncludeCurrent())
tf.push_back(ExtraRotamers(0, 1, 1))
tf.push_back(ExtraRotamers(0, 2, 1))

prevent = PreventRepackingRLT() # No repack, no design
repack = RestrictToRepackingRLT() # No design
tf.push_back(OperateOnResidueSubset(prevent, static_res))
tf.push_back(OperateOnResidueSubset(repack, packable_res))
tf.push_back(OperateOnResidueSubset(repack, ligand))

fd = FastDesign()
fd.set_scorefxn(sf)
fd.set_task_factory(tf)

filnam = 'start_models/{}.pdb'.format(args.ligand)
pose = pose_from_pdb(filnam)

jnam = 'fixed_designs/{}_designed_{}'.format(args.ligand, args.job)
jd = PyJobDistributor(jnam, 1, sf)
while not jd.job_complete:
	pp = Pose(pose)
	fd.apply(pp)
	jd.output_decoy(pp)

