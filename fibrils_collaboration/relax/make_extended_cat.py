from pyrosetta import *
init()
pose=pose_from_pdb('fibrils_collaboration/relax/20190201/htra1_protease_relaxed_2_6.pdb.gz')
sp = pose.split_by_chain()[2]
pep=pose_from_sequence('GMFKLVVAKPPS') # P6 and P' from Ehrmann WT/soluble AA dist most frequent
for i in range(1,sp.total_residue()+1):
    pep.set_phi(i+1, sp.phi(i))
    pep.set_psi(i+1, sp.psi(i))
    pep.set_chi(1, i+1, sp.chi(1,i))
pep.set_phi(2,-180)
pep.set_chi(2,2,sp.chi(2,1))
pep.set_chi(3,2,sp.chi(3,1))
pep.set_chi(2,3,sp.chi(2,2))
pep.set_chi(2,4,sp.chi(2,3))
pep.set_chi(3,4,sp.chi(3,3))
pep.set_chi(4,4,sp.chi(4,3))
pep.set_chi(2,5,sp.chi(2,4))
pep.dump_pdb('fibrils_collaboration/extended_pep.pdb')
"""
load into pymol
alter chain to B
pair fit to align
manually adjust dihedrals of downstream chain to align near the b-sheet (res 200)
add in constraints block, adjusting for renumbered chain B
"""