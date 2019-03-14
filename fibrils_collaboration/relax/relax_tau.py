#!/usr/bin/python
from pyrosetta import *
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover

init('-relax:thorough')

sf = SymmetricScoreFunction()
sf.add_weights_from_file('ref2015_cst')

fr = FastRelax(sf, "NO CST RAMPING")
fr.set_scorefxn(sf)
fr.constrain_coords(True)
fr.constrain_relax_to_start_coords(True)
fr.ramp_down_constraints(False)
fr.coord_constrain_sidechains(True)

sfsm = SetupForSymmetryMover('tau_fibril.sym')

pose = pose_from_pdb('5o3l_INPUT.pdb')
sfsm.apply(pose)

jd = PyJobDistributor('sym_tau_relax/tau_relax', 50, sf)
while not jd.job_complete:
	p = Pose()
	p.assign(pose)
	fr.apply(p)
	jd.output_decoy(p)