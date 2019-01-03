#!/usr/bin/python
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts

init('-enzdes::cstfile htra1_cat_general.cst -cst_fa_weight 1.0 -run:preserve_header')

sf = create_score_function('ref2015_cst')

fr = FastRelax()
fr.set_scorefxn(sf)
fr.constrain_coords(True)
fr.constrain_relax_to_start_coords(True)
fr.coord_constrain_sidechains(True)

cstm = AddOrRemoveMatchCsts()
cstm.set_cst_action(ADD_NEW)

pose = pose_from_pdb('cat_pre_relax.pdb')
cstm.apply(pose)

jd = PyJobDistributor('cat_relax/cat_relax', 5, get_fa_scorefxn())
while not jd.job_complete:
	pp = Pose()
	pp.assign(pose)
	fr.apply(pp)
	jd.output_decoy(pp)