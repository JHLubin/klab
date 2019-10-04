from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import ExtraRotamers, IncludeCurrent, OperateOnResidueSubset, PreventRepackingRLT, RestrictToRepackingRLT, RestrictAbsentCanonicalAASRLT
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector, ChainSelector, NeighborhoodResidueSelector, NotResidueSelector, OrResidueSelector, ResidueIndexSelector
init('-run:preserve_header')
sf=get_fa_scorefxn()
pose = pose_from_pdb('pdz_relaxed_2.pdb')
for i in [382, 386, 388, 389, 394, 415, 416, 444, 446]:
	print(i, pose.pdb_info().pdb2pose('A',i))
"""
382 7
386 11
388 13
389 14
394 19
415 40
416 41
444 69
446 71
"""	
res_changes = {
'wt': {106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'Y382W': {7:'W', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'R386V': {11:'V', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'R386V-M388R-K394V': {11:'V', 13:'R', 19:'V', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'I415K': {40:'K', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'I415K-E416Y': {40:'K', 41:'Y', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'R386V-M388R-K394V-T415K-E416Y': {11:'V', 13:'R', 19:'V', 40:'K', 41:'Y', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'Y382W-R386V-M388R-K394V-T415K-E416Y': {7:'W', 11:'V', 13:'R', 19:'V', 40:'K', 41:'Y', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'S389D-S444R-N446Q': {14:'D', 69:'R', 71:'Q', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'},
'N446E': {71:'E', 106:'Q', 107:'D', 108:'Y', 109:'E', 110:'P', 111:'E', 112:'A'}
}
start_pdb = {
'Y382W': 'pdz_relaxed_1.pdb_designed_62.pdb',
'R386V': 'pdz_relaxed_2.pdb_designed_13.pdb',
'R386V-M388R-K394V': 'pdz_relaxed_2.pdb_designed_13.pdb',
'I415K': 'pdz_relaxed_2.pdb_designed_13.pdb',
'I415K-E416Y': 'pdz_relaxed_2.pdb_designed_13.pdb',
'R386V-M388R-K394V-T415K-E416Y': 'pdz_relaxed_2.pdb_designed_13.pdb',
'Y382W-R386V-M388R-K394V-T415K-E416Y': 'pdz_relaxed_2.pdb_designed_13.pdb',
'S389D-S444R-N446Q': 'pdz_relaxed_3.pdb_designed_47.pdb',
'N446E': 'pdz_relaxed_2.pdb_designed_9.pdb'
}
for des, muts in res_changes.items():
	pp = pose_from_pdb('pdz_designs/examples/' + start_pdb[des])
	peptide = ChainSelector('B')
	shell = NeighborhoodResidueSelector()
	shell.set_focus_selector(peptide)
	shell.set_include_focus_in_subset(True)
	shell.set_distance(12)
	mobile_residues = OrResidueSelector()
	tf = TaskFactory()
	tf.push_back(IncludeCurrent())
	tf.push_back(ExtraRotamers(0, 1, 1))
	tf.push_back(ExtraRotamers(0, 2, 1))
	for r, aa in muts.items():
		res_selection = ResidueIndexSelector(str(r))
		restriction = RestrictAbsentCanonicalAASRLT()
		restriction.aas_to_keep(aa.upper())
		tf.push_back(OperateOnResidueSubset(restriction,res_selection))
		mobile_residues.add_residue_selector(res_selection)
	packable = AndResidueSelector(shell, NotResidueSelector(mobile_residues))
	tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), packable))
	not_changing = NotResidueSelector(shell)
	tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), not_changing))
	pt=tf.create_task_and_apply_taskoperations(pose)
	prm=PackRotamersMover(sf, pt)
	prm.apply(pp)
	pp.dump_pdb('pdz_designs/design_models/pdz_' + des + '.pdb')


	# Run analyze_design_decoys.fix_file()
	# python relax_new_pdb.py $i -od fibrils_collaboration/pdz_designs/design_models/round_2/ -cst fibrils_collaboration/htra1_pdz.cst -ccw 0 -n 10