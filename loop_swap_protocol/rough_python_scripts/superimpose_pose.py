from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import PerResidueRMSDMetric
from pyrosetta.rosetta.core.scoring import rmsd_atoms
from pyrosetta.rosetta.core.scoring import superimpose_pose
init()
pose=pose_from_pdb('loop_swap_protocol/tev.pdb')
ky9 = pose_from_pdb('../dali_data/pdb/ky/pdb1ky9.ent.gz')
pmm = PyMOLMover()
pmm.apply(pose)
pmm.apply(ky9)
qris = ResidueIndexSelector(str(pose.pdb_info().pdb2pose('A',142))+","+str(pose.pdb_info().pdb2pose('A',154)))
sris = ResidueIndexSelector(str(ky9.pdb_info().pdb2pose('B',201))+","+str(ky9.pdb_info().pdb2pose('B',213)))
prmsd = PerResidueRMSDMetric()
prmsd.set_rmsd_type(rmsd_atoms.rmsd_protein_bb_ca)
prmsd.set_residue_selector_reference(qris)
prmsd.set_residue_selector(sris)
prmsd.set_comparison_pose(pose)
amap = prmsd.create_atom_id_map(ky9)
superimpose_pose(ky9, pose, amap)
pmm.apply(ky9)

			# Make residue selectors for loop termini with overlap regions
			stemp = '{}-{},{}-{}'
			qs = stemp.format(qn - self.n_overlap_size, qn, 
								qc, qc + self.c_overlap_size)
			ss = stemp.format(sn - self.n_overlap_size, sn, 
								sc, sc + self.c_overlap_size)
			query_selector = ResidueIndexSelector(qs)
			subject_selector = ResidueIndexSelector(ss)

rmsd = RMSDMetric()


prmsd.set_residue_selector_reference?
psubset = Pose(pose(1,50))
Pose?
ky9=pose_from_pdb('loop_swap_protocol/aligned_pdbs/1KY9.pdb')
tev_inf = pose.pdb_info()
ky9_inf=ky9.pdb_info()
qloop = Pose(pose, tev_inf.pdb2pose('A',142),tev_inf.pdb2pose('A',154))
sloop = Pose(ky9, ky9_inf.pdb2pose('A',201),ky9_inf.pdb2pose('A',213))
pmm = PyMOLMover()
pmm.apply(qloop)
pmm.apply(sloop)
pmm.apply(qloop)
pmm.apply(sloop)
pmm.apply(pose)
rmsd.set_comparison_pose?
rmsd.set_comparison_pose(qloop)
rmsd.calculate(sloop)
rmsd.set_rmsd_type?
rmsd.set_rmsd_type('CA')
rmsd.set_rmsd_type('rmsd_protein_bb_ca')
rmsd_atoms?
rmsd.set_rmsd_type(rmsd_atoms.rmsd_protein_bb_ca)
rmsd.calculate(sloop)
prmsd.set_comparison_pose(qloop)
prmsd.calculate(sloop)
prmsd.set_rmsd_type(rmsd_atoms.rmsd_protein_bb_ca)
prmsd.calculate(sloop)
sloop.total_residue()
print(sloop)
print(qloop)
ris = ResidueIndexSelector('1,13')
ris = ResidueIndexSelector(str(ky9_inf.pdb2pose('A',201))+","+str(ky9_inf.pdb2pose('A',213)))
prmsd.set_residue_selector_reference(ris)
sris = ResidueIndexSelector(str(ky9_inf.pdb2pose('A',201))+","+str(ky9_inf.pdb2pose('A',213)))
qris = ResidueIndexSelector(str(tev_inf.pdb2pose('A',142))+","+str(tev_inf.pdb2pose('A',154)))
prmsd.set_residue_selector_reference(qris)
prmsd.set_residue_selector(sris)
prmsd.set_comparison_pose(pose)
amap = prmsd.create_atom_id_map(ky9)
superimpose_pose(ky9, pose, amap)
pmm.apply(ky9)
superimpose_pose(ky9, *pose, amap)
ky9_b = pose_from_pdb('../dali_data/pdb/ky/pdb1ky9.ent.gz')
print(ky9_b)
sris_b = ResidueIndexSelector(str(ky9_b.pdb_info().pdb2pose('B',201))+","+str(ky9_b.pdb_info().pdb2pose('B',213)))
prmsd.set_residue_selector(sris_b)
amap = prmsd.create_atom_id_map(ky9_b)
pmm.apply(ky9_b)
superimpose_pose(ky9_b, pose, amap)
pmm.apply(ky9_b)
ky9_b = pose_from_pdb('../dali_data/pdb/ky/pdb1ky9.ent.gz')