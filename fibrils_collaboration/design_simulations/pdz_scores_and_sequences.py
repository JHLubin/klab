from pyrosetta import *
from pyrosetta.rosetta.core.simple_metrics.metrics import SequenceMetric
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector, ChainSelector, InterGroupInterfaceByVectorSelector, NotResidueSelector, OrResidueSelector, ResidueIndexSelector
from glob import glob

init('-mute all')
pose=pose_from_pdb('../../cat_relaxed.pdb')
pdbs = glob('*designed*.pdb')
pdbs.sort()
cha = ChainSelector('A')
a = SequenceMetric(cha)
orig_seq = a.calculate(pose)

for pdb in pdbs:
p = pose_from_pdb(pdb)
p_seq = a.calculate(p)
muts = []
for n, i in enumerate(orig_seq):
if p_seq[n] != i:
muts.append(i+p.pdb_info().pose2pdb(n+1).split()[0]+p_seq[n])

print(pdb, ' '.join(muts))

pose = pose_from_pdb('../start_proteases/HCV.pdb')
chb = ChainSelector('B')
b = SequenceMetric(chb)
pep_seq = b.calculate(pose)
for pdb in pdbs:
p = pose_from_pdb(pdb)
print(pdb, b.calculate(p)[1:6])