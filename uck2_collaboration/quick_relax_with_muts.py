import argparse
from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, NotResidueSelector, OrResidueSelector
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	RestrictToRepackingRLT, RestrictToRepacking, RestrictAbsentCanonicalAASRLT, OperateOnResidueSubset
from pyrosetta.rosetta.protocols.relax import FastRelax

parser = argparse.ArgumentParser()
parser.add_argument('--ligand', help='Which ligand?')
parser.add_argument('--mutant', nargs=2, action='append', help='site residue (res 65 is 39 in pose)')
parser.add_argument('--n', type=int, help='How many decoys?')
args = parser.parse_args()

opts = '-use_input_sc -ex1 -ex2 -ignore_zero_occupancy false -nblist_autoupdate  -extra_res_fa params/{}.params'.format(args.ligand)
init(opts)

sf=create_score_function('ref2015_cst')
fr=FastRelax()
fr.set_scorefxn(sf)
fr.constrain_relax_to_start_coords(True)

tf=TaskFactory()
if args.mutant:
	for i in args.mutant:
		res_selection = ResidueIndexSelector(str(i[0]))
		restriction = RestrictAbsentCanonicalAASRLT()
		restriction.aas_to_keep(i[1].upper())
		tf.push_back(OperateOnResidueSubset(restriction, res_selection))
		tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(),NotResidueSelector(res_selection)))
else:
	tf.push_back(RestrictToRepacking)

fr.set_task_factory(tf)

filnam = 'start_models/{}.pdb'.format(args.ligand)
pose = pose_from_pdb(filnam)

jnam = 'check_muts/{}_{}{}'.format(args.ligand, *args.mutant[0])
jd = PyJobDistributor(jnam, args.n, sf)
while not jd.job_complete:
	pp = Pose(pose)
	fr.apply(pp)
	jd.output_decoy(pp)