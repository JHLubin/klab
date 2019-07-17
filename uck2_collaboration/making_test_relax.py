cdkl; cd uck2_collaboration/; ipython
lig = '5cpu'
from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, NotResidueSelector
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import OperateOnResidueSubset, PreventRepackingRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.protocols.relax import FastRelax
#opts = '-use_input_sc -ex1 -ex2 -ignore_zero_occupancy false -nblist_autoupdate  -extra_res_fa params/{}.params'.format(lig)
opts = '-use_input_sc -ignore_zero_occupancy false -nblist_autoupdate  -extra_res_fa params/{}.params'.format(lig)
init(opts)
sf=get_fa_scorefxn()
fr=FastRelax()
fr.set_scorefxn(sf)
close_res = ResidueIndexSelector('34,35,36,37,38,39,40,41,55,56,57,58,59,60,61,62,63,64,65,68,86,87,88,89,90,91,93,108,109,110,111,112,113,148,205')
tf = TaskFactory()
prevent = PreventRepackingRLT()
repack = RestrictToRepackingRLT()
tf.push_back(OperateOnResidueSubset(repack, close_res))
tf.push_back(OperateOnResidueSubset(prevent, NotResidueSelector(close_res)))
fr.set_task_factory(tf)
pose=pose_from_pdb('start_models/{}.pdb'.format(lig))
jd=PyJobDistributor('test_relax/{}'.format(lig),3,sf)
while not jd.job_complete:
pp=Pose(pose)
fr.apply(pp)
jd.output_decoy(pp)

