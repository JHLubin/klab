from pyrosetta import *
init()
pose=pose_from_pdb('orig_cat_relaxed.pdb')
sp = pose.split_by_chain()[2]
pep=pose_from_sequence('TMFKLVTTTTTT')
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
pep.dump_pdb('extended_pep.pdb')
"""
load into pymol
alter chain to B
pair fit to align
manually adjust dihedrals of downstream chain to align near the b-sheet (res 200)
"""