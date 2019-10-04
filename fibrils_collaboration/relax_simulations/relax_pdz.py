#!/usr/bin/python
from glob import glob
from os.path import basename, join
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax

init()

sf = create_score_function('ref2015_cst')

fr = FastRelax()
fr.set_scorefxn(sf)
fr.constrain_coords(True)
fr.constrain_relax_to_start_coords(True)
fr.coord_constrain_sidechains(True)

structs = glob('pdz_nmr/*.pdb')
structs.sort()

for s in structs:
	pose = pose_from_pdb(s)

	jdname = join('pdz_relax/', basename(s).rstrip('.pdb'))

	jd = PyJobDistributor(jdname, 20, sf)
	while not jd.job_complete:
		pp = Pose()
		pp.assign(pose)
		fr.apply(pp)
		jd.output_decoy(pp)