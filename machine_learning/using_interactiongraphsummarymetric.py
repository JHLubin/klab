from pyrosetta import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import IncludeCurrent, OperateOnResidueSubset, PreventRepackingRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, NeighborhoodResidueSelector, NotResidueSelector, OrResidueSelector
from pyrosetta.rosetta.protocols.quantum_annealing import InteractionGraphSummaryMetric
init()

# Starting pose
pose=pose_from_pdb('HCV.pdb')

# Score function
sf = get_fa_scorefxn()

#Selections
peptide = ChainSelector('B')

shell_1 = NeighborhoodResidueSelector()
shell_1.set_distance(8)
shell_1.set_focus_selector(peptide)
shell_1.set_include_focus_in_subset(False)

shell_2 = NeighborhoodResidueSelector()
shell_2.set_distance(4)
shell_2.set_focus_selector(OrResidueSelector(peptide, shell1))
shell_2.set_focus_selector(OrResidueSelector(peptide, shell_1))
shell_2.set_include_focus_in_subset(False)

mobile = OrResidueSelector()
mobile.add_residue_selector(peptide)
mobile.add_residue_selector(shell_1)
mobile.add_residue_selector(shell_2)

immobile = NotResidueSelector(mobile)

# Making packer task
tf = TaskFactory()
prevent = PreventRepackingRLT()
repack = RestrictToRepackingRLT()
tf.push_back(IncludeCurrent())
tf.push_back(OperateOnResidueSubset(repack, peptide))
tf.push_back(OperateOnResidueSubset(repack, shell_2))
tf.push_back(OperateOnResidueSubset(prevent, immobile))

# Interaction graph summary
igsm = InteractionGraphSummaryMetric()
igsm.set_task_factory(tf)
igsm.set_scorefunction(sf)
igsm.set_short_version(True)

# Output file
with open('interaction_graph_summary.txt', 'w') as o:
	o.write(igsm.calculate(pose))