cdkl
cd fibrils_collaboration
ipython
pdb_name = 'pdz_designs/design_models/pdz_Y382W-R386V-M388R-K394V-T415K-E416Y.pdb'
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
ros_opts = '-ex1 -ex2  -use_input_sc -flip_HNQ -enzdes::cstfile {} -run:preserve_header'
init(ros_opts.format('htra1_pdz.cst'))
sf = get_fa_scorefxn()
sf.set_weight(ScoreType.atom_pair_constraint, 1)
sf.set_weight(ScoreType.coordinate_constraint, 1)
sf.set_weight(ScoreType.angle_constraint, 1)
sf.set_weight(ScoreType.dihedral_constraint, 1)
fr = FastRelax()
fr.set_scorefxn(sf)
pose = pose_from_pdb(pdb_name)
cstm = AddOrRemoveMatchCsts()
cstm.set_cst_action(ADD_NEW)
cstm.apply(pose)
jdname = pdb_name.replace('.pdb', '_relaxed')
jd = PyJobDistributor(jdname, 10, sf)
while not jd.job_complete:
pp = Pose(pose)
fr.apply(pp)
jd.output_decoy(pp)

