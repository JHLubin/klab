from pyrosetta import *
from pyrosetta.rosetta.protocols.grafting.simple_movers import DeleteRegionMover, KeepRegionMover, ReplaceRegionMover
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.grafting import CCDEndsGraftMover
init()
pose=pose_from_pdb('a_to_s_ly104_WT.pdb')
p1=Pose()
p1.assign(pose)
drm=DeleteRegionMover()
ri=ResidueIndexSelector('135-138')
drm.set_residue_selector(ri)
drm.apply(p1)
p2=Pose()
p2.assign(pose)
krm=KeepRegionMover()
krm.start(135)
krm.end(138)
krm.apply(p2)
p3=Pose()
p3.assign(p1)
p3.append_polymer_residue_after_seqpos(p2.residue(1),134,False)
p3.append_polymer_residue_after_seqpos(p2.residue(2),135,False)
p3.append_polymer_residue_after_seqpos(p2.residue(3),136,False)
p3.append_polymer_residue_after_seqpos(p2.residue(4),137,False)
p3.append_polymer_residue_after_seqpos(p2.residue(5),138,False)
p3.assign(pose)
ri=ResidueIndexSelector('133-139')
drm.set_residue_selector(ri)
drm.apply(p3)
hrv=pose_from_pdb('serine_proteases/combined_list/2hrv.pdb')
loop=Pose()
loop.assign(hrv)
krm.start(72)
krm.end(97)
krm.apply(loop)
pm.apply(p3)
pm.apply(lop)
pm.apply(loop)
begin=132
for i in range(1,loop.total_residue()+1):
    p3.append_polymer_residue_after_seqpos(loop.residue(i),begin,False)
    begin+=1
rrm=ReplaceRegionMover(hrv, 72, 1, 26, False)
aaaa=pose_from_sequence('A'*100)
rrm.apply(aaaa)
hrv=pose_from_pdb('serine_proteases/combined_list/2hrv.pdb')
krm.start(74); krm.end(95)
krm.apply(hrv)
pose=pose_from_pdb('a_to_s_ly104_WT.pdb')
ccdgm=CCDEndsGraftMover(135,138,hrv,2,2)
ccdgm.apply(pose)

poses={}
anchors=[[134,139], [135,138], [136,137]]
loops=[[72,97], [73,96], [74,95], [75,94], [76,93]]
overlaps=range(7)
for a in anchors:
	for l in loops:
		for o in overlaps:
			name=str(a) + str(l) + str(o)
			name=name.replace('[','').replace(']','_').replace(',','-')
			name=name.replace(' ', '')
			name='calibrate_CCDEndsGraftMover/'+name+'.pdb'
			if name not in made:
				try:
					hrv=pose_from_pdb('2hrv.pdb')
					krm.start(l[0]); krm.end(l[1])
					krm.apply(hrv)
					pose=pose_from_pdb('a_to_s_ly104_WT.pdb')
					ccdgm=CCDEndsGraftMover(a[0],a[1],hrv,o,o)
					ccdgm.apply(pose)
					pose.dump_pdb('calibrate_CCDEndsGraftMover/'+name+'.pdb')