import loop_swap_protocol.protein_loop_alignment as pla
from math import atan, pi, radians, sqrt
from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.scoring.constraints import \
    AngleConstraint, AtomPairConstraint, DihedralConstraint
from pyrosetta.rosetta.core.scoring.func import \
    CircularHarmonicFunc, HarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
init()
tev=pose_from_pdb('loop_swap_protocol/tev.pdb')
tev.pdb_info().pdb2pose('A',142)
tev.pdb_info().pdb2pose('A',154)
qris = ResidueIndexSelector('134-135,147-148')
ky9 = pose_from_pdb('../dali_data/pdb/ky/pdb1ky9.ent.gz')
ky9.pdb_info().pdb2pose('B',201)
ky9.pdb_info().pdb2pose('B',213)
sris = ResidueIndexSelector('466-467,479-480')
pla.align_protein_sections(tev, qris, ky9, sris)
loop = Pose(ky9, 466,480)
loop.append_pose_by_jump(pose_from_sequence('G'),1)
size_l = loop.total_residue()
n_ca = AtomID(loop.residue(2).atom_index('CA'), 2)
c_ca = AtomID(loop.residue(size_l - 2).atom_index('CA'), size_l -2)
x_ca = AtomID(loop.residue(size_l).atom_index('CA'), size_l)
x_n = AtomID(loop.residue(size_l).atom_index('N'), size_l)
x_c = AtomID(loop.residue(size_l).atom_index('C'), size_l)
r_ca = AtomID(loop.residue(1).atom_index('CA'), 1)
int_ca_dist = pla.get_distance(pla.find_res_ca_coords(loop,1),pla.find_res_ca_coords(loop,size_l-1)) / 2
out_dist1 = sqrt(int_ca_dist**2 + 20**2)
out_dist2 = sqrt(int_ca_dist**2 + 21.5**2)
dist_hf1 = HarmonicFunc(out_dist1, 0.1)
dist_hf2 = HarmonicFunc(out_dist2, 0.1)
dcst1 = AtomPairConstraint(n_ca, x_n, dist_hf1)
dcst2 = AtomPairConstraint(c_ca, x_n, dist_hf1)
dcst3 = AtomPairConstraint(n_ca, x_ca, dist_hf2)
dcst4 = AtomPairConstraint(c_ca, x_ca, dist_hf2)
out_ang = atan(20/int_ca_dist)
ang_hf1 = HarmonicFunc(out_ang,1)
ang_hf2 = HarmonicFunc(pi/2 + out_ang, 1)
acst1 = AngleConstraint(c_ca, n_ca, x_n, ang_hf)
#acst2 = AngleConstraint(n_ca, c_ca, x_n, ang_hf)
acst2 = AngleConstraint(n_ca, x_n, x_ca, ang_hf2)
chf1 = CircularHarmonicFunc(0, 0.01)
chf2 = CircularHarmonicFunc(pi, 0.01)
chf3 = CircularHarmonicFunc(pi/2, 0.01)
ccst1 = DihedralConstraint(r_ca, c_ca, n_ca, x_n, chf1)
ccst2 = DihedralConstraint(c_ca, n_ca, x_n, x_ca, chf2)
ccst3 = DihedralConstraint(n_ca, x_n, x_ca, x_c, chf3)
#ccst4 = DihedralConstraint(n_ca, x_c, x_ca, x_n, chf3)
#constratint_set = [dcst1, dcst2, acst1, acst2, ccst1, ccst2, ccst3, ccst4]
constratint_set = [dcst, acst1, acst2, ccst1, ccst2, ccst3]
for c in constratint_set:
    loop.add_constraint(c)
sf = create_score_function('ref2015_cst')
mm=MoveMap()
mm.set_jump(True)
mm.set_bb(False)
mm.set_chi(False)
minmov = MinMover()
minmov.score_function(sf)
minmov.min_type('lbfgs_armijo_nonmonotone')
minmov.movemap(mm)

In [262]: int_ca_dist
Out[262]: 5.188580947619549

In [263]: out_dist
Out[263]: 20.662075700422758

In [264]: degrees(out_ang)
Out[264]: 75.4564298506566

In [265]: degrees(out_ang+pi/2)
Out[265]: 165.45642985065658

In [266]: 180-degrees(out_ang+pi/2)
Out[266]: 14.543570149343424

yol = pose_from_pdb('../dali_data/pdb/yo/pdb2yol.ent.gz')
yol.pdb_info().pdb2pose('A',1028)
yol.pdb_info().pdb2pose('A',1033)
tev.pdb_info().pdb2pose('A',25)
tev.pdb_info().pdb2pose('A',28)
qris = ResidueIndexSelector('24-25,28-29')
sris = ResidueIndexSelector('70-71,76-77')
pla.align_protein_sections(tev, qris, yol, sris)
loop = Pose(yol, 70,77)
int_ca_dist = pla.get_distance(pla.find_res_ca_coords(loop,2),pla.find_res_ca_coords(loop,size_l-2)) / 2
out_dist1 = sqrt(int_ca_dist**2 + 20**2)
out_dist2 = sqrt(int_ca_dist**2 + 21.5**2)
out_ang = atan(20/int_ca_dist)
dist_hf1 = HarmonicFunc(out_dist1, 0.5)
dist_hf2 = HarmonicFunc(out_dist2, 0.5)
dcst1 = AtomPairConstraint(n_ca, x_n, dist_hf1)
dcst2 = AtomPairConstraint(c_ca, x_n, dist_hf1)
ang_hf1 = HarmonicFunc(out_ang,0.1)
ang_hf2 = HarmonicFunc(pi/2 + out_ang, 0.5)
acst1 = AngleConstraint(n_ca, x_n, x_ca, ang_hf2)
acst2 = AngleConstraint(c_ca, x_n, x_ca, ang_hf2)
chf1 = CircularHarmonicFunc(0, 0.01)
chf2 = CircularHarmonicFunc(pi, 0.01)
chf3 = CircularHarmonicFunc(pi/2, 0.01)
ccst1 = DihedralConstraint(r_ca, c_ca, n_ca, x_n, chf1)
ccst2 = DihedralConstraint(c_ca, n_ca, x_n, x_ca, chf1)
ccst3 = DihedralConstraint(n_ca, x_n, x_ca, x_c, chf3)
pp=Pose(loop)
for c in [dcst1, dcst2, acst2, ccst1, ccst2, ccst3]:
    pp.add_constraint(c)
minmov.apply(pp); pmm.apply(pp); sf.show(pp)