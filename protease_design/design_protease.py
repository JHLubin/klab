#!/usr/bin/python
""" 
PyRosetta4, Python 2.7
Joseph Lubin, 2017
"""
import argparse
from math import sqrt
from os import makedirs
from os.path import basename, isdir, join
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.select.residue_selector import \
	AndResidueSelector, ChainSelector, InterGroupInterfaceByVectorSelector,\
	NotResidueSelector, OrResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.core.pack.task import parse_resfile, TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	OperateOnResidueSubset, PreventRepackingRLT, \
	RestrictAbsentCanonicalAASRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.numeric import xyzVector_double_t
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.teaching import MinMover, PackRotamersMover, SimpleThreadingMover
from sys import exit

def parse_args():
	info = "Design a protease around a peptide sequence"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-s", "--start_struct", required=True,
		default='a_to_s_ly104_WT.pdb', help="Pick starting PDB")
	parser.add_argument("-od", "--out_dir", required=True,
		help="Name an output directory for decoys")
	parser.add_argument("-cseq", "--cut_peptide_sequence", type=str, 
		action='append', default=[], help="List cleaved peptide sequences or \
		 provide a list file")
	parser.add_argument("-useq", "--uncut_peptide_sequence", type=str, 
		action='append', default=[], help="List uncleaved peptide sequences \
		or provide a list file")
	parser.add_argument("-cr", "--cat_res", type=int, nargs='+', 
		default=[72, 96, 154], help="The catalytic residues of the protease, \
		excluded from design (defaults are 72, 96, and 154, for HCV)")
	parser.add_argument("-cons", "--constraints", type=str, 
		default='ly104.cst', help="Pick constraints file")
	parser.add_argument("-des_pep", "--design_peptide", action="store_true", 
		help="Option to allow design on the regognition region of the peptide")
	parser.add_argument("-th", "--thread", action="store_true", default=False,
		help="Option to create threaded models. Use the first time running.")
	parser.add_argument("-d", "--number_decoys", type=int, default=10, 
		help="How many decoys should be made? (Default is 10.)")
	parser.add_argument("-m", "--method", type=str, default='fastdesign', 
		choices=['fastdesign', 'custom'],
		help="What design protocol do you want to use?")
	parser.add_argument("-trmf", "--target_res_mutation_file", type=str, 
		default=None, help="File of format: 'res_number    design_options' \
		to specify further design restrictions. Useful when picking between \
		identified favorable mutations.")
	args = parser.parse_args()
	return args


def init_opts(cst_file='ly104.cst'):
	""" Produces a list of init options for PyRosetta, including cst file """
	ros_opts = '-mute all -enzdes::cstfile ' + cst_file
	ros_opts += ' -cst_fa_weight 1.0 -run:preserve_header -out:pdb_gz'
	return ros_opts


def readfile(file_name):
	""" Opens a file in read-mode and returns a list of the text lines """
	with open(file_name, 'r') as r:
		lines = r.readlines()

	return lines


######### Threading and getting sequences ####################################
def get_seq_list(seq_arg):
	"""
	Takes an argument that can include individual peptide sequences or file(s)
	containing a list of sequences, and returns a list of sequences. 
	Distinguishes between files and sequences by the presence of a dot (.).
	"""
	pep_sequences = []
	for inp in seq_arg:
		if '.' in inp:
			# If input is a file
			lis = readfile(inp) 
			if len(lis) == 1:
				# If all sequences are listed horizontally on one line
				# rather than one per line, rearrange
				lis = lis[0].split()

			for i in lis:
				pep_sequences.append(i.strip())

		else:
			# Sequence was typed directly into the argument
			pep_sequences.append(inp.strip())

		for n, s in enumerate(pep_sequences):
			if len(s) == 5: # Hard code for HCV
				pep_sequences[n] = 'A' + s + 'CSMHL'
			elif len(s) == 6: # Hard code for HCV
				pep_sequences[n] = 'A' + s + 'SMHL'
			else:
				assert len(s) == 11

	return pep_sequences


def quick_thread(destination, pose, sequences, cleaved=False, make=False):
	""" 
	Threads a set of sequences onto the peptide portion of the given PDB file,
	outputting a threaded PDB file for each sequence.
	Function is presently hard-coded for this application.
	"""
	thread_files = []
	for seq in sequences:
		# Naming model
		pdbname = 'cleaved_ly104_wt_' + seq + '.pdb.gz'
		if not cleaved:
			pdbname = 'un' + pdbname

		out_name = join(destination, pdbname)
		thread_files.append(out_name)

		if make:
			print(out_name)
			# Threading peptide sequences
			tm = SimpleThreadingMover(seq, 197)
			threaded_pose = Pose()
			threaded_pose.assign(pose)
			tm.apply(threaded_pose)
			threaded_pose.dump_pdb(out_name)

	return thread_files


######### Residue selection ##################################################
def mutable_residues_selector(design_peptide=False, single_pep_res=None):
	"""
	Selects the residues in a shell around the peptide using the 
	InterGroupInterfaceByVectorSelector residue selector
	Presently hard coded for HCV protease.
	Decided to design just around peptide, not pep + cat, and only the
	six mutable residues in the peptide.
	"""
	# Selecting regions
	protease = ChainSelector("A")
	catalytic = ResidueIndexSelector('72,96,154')
	if single_pep_res:
		variable_pep_res = ResidueIndexSelector(str(single_pep_res))
	else:
		variable_pep_res = ResidueIndexSelector("198-203")

	# Making positive residue selector
	rs = InterGroupInterfaceByVectorSelector()
	rs.group1_selector(protease) # Protease
	rs.group2_selector(variable_pep_res) # Peptide recognition region
	rs.nearby_atom_cut(6)
	rs.vector_dist_cut(8)

	# Excluding the catalytic residues
	limit_selection = AndResidueSelector()
	limit_selection.add_residue_selector(NotResidueSelector(catalytic))
	limit_selection.add_residue_selector(rs)

	# If the peptide sequence is mutable
	if design_peptide:
		return limit_selection

	# If only the protease is designable
	else: 
		# Setting up exclusion of catalytic and peptide residues
		exclusive_selection = AndResidueSelector()
		exclusive_selection.add_residue_selector(limit_selection)
		exclusive_selection.add_residue_selector(ChainSelector("A"))
		return exclusive_selection


def packable_residues_selector(mutable_selector):
	"""
	Selects the shell of neighbor residues to repack
	Presently hard coded for HCV protease.
	"""
	# Selecting regions
	mutable = mutable_selector
	not_mutable = NotResidueSelector(mutable_selector)
	catalytic = ResidueIndexSelector('72,96,154')
	peptide = ChainSelector('B')

	pep_not_mutable = AndResidueSelector()
	pep_not_mutable.add_residue_selector(peptide)
	pep_not_mutable.add_residue_selector(not_mutable)

	# Selecting residues near mutable shell
	near_mutable = InterGroupInterfaceByVectorSelector()
	near_mutable.group1_selector(not_mutable)
	near_mutable.group2_selector(mutable)
	near_mutable.nearby_atom_cut(4)
	near_mutable.vector_dist_cut(4)

	# Selecting residues near the peptide
	near_pep = InterGroupInterfaceByVectorSelector()
	near_pep.group1_selector(not_mutable) # Protease
	near_pep.group2_selector(peptide) # Peptide recognition region
	near_pep.nearby_atom_cut(10)
	near_pep.vector_dist_cut(12)

	# Combining selections for near peptide and near mutable residues
	wide_set = OrResidueSelector()
	wide_set.add_residue_selector(near_mutable)
	wide_set.add_residue_selector(near_pep)	

	# Setting up exclusion of catalytic and mutable residues
	limit_selection = AndResidueSelector()
	limit_selection.add_residue_selector(NotResidueSelector(catalytic))
	limit_selection.add_residue_selector(wide_set)
	limit_selection.add_residue_selector(not_mutable)

	# Add back in the peptide
	expand_selection = OrResidueSelector()
	expand_selection.add_residue_selector(limit_selection)
	expand_selection.add_residue_selector(pep_not_mutable)

	return expand_selection


def other_residues_selector(mutable_selector, packable_selector):
	""" Selects the residues that are not designable or repackable """
	all_mobile_res = OrResidueSelector()
	all_mobile_res.add_residue_selector(mutable_selector)
	all_mobile_res.add_residue_selector(packable_selector)

	other_res_selector = NotResidueSelector(all_mobile_res)

	return other_res_selector


def movemap_packable_residues_selector(peptides, mutables, packables):
	""" 
	Combines selectors for all residues for which a MoveMap should 
	set chi true.
	"""
	combined = OrResidueSelector()
	combined.add_residue_selector(peptides)
	combined.add_residue_selector(mutables)
	combined.add_residue_selector(packables)

	return combined


def select_residues(design_peptide=False):
	""" Makes residue selectors for HCV protease sections """
	residue_selectors = {}

	catalytic = ResidueIndexSelector('72,96,154')
	residue_selectors['catalytic'] = catalytic

	protease = ChainSelector("A")
	residue_selectors['protease'] = protease

	peptide = ChainSelector("B")
	residue_selectors['peptide'] = peptide

	mutable = mutable_residues_selector(design_peptide)
	residue_selectors['mutable'] = mutable

	packable = packable_residues_selector(mutable)
	residue_selectors['packable'] = packable
	
	other = other_residues_selector(mutable, packable)
	residue_selectors['other'] = other

	mm_pack = movemap_packable_residues_selector(peptide, mutable, packable)
	residue_selectors['movemap_pack'] = mm_pack

	return residue_selectors


def selector_to_list(pose, selector):
	""" Converts a selector output vector to a list of selected residues """
	selection_vector = selector.apply(pose)
	selection_list = []
	for i in range(len(selection_vector)): 
		if selection_vector[i+1]==1:
			selection_list.append(i+1)

	return selection_list 


######### Design #############################################################
def apply_constraints(pose):
	""" Applies the constraints form the input CST file to a pose """
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)
	return pose


def make_fold_tree():
	"""
	Make a fold tree that connects the first catalytic residue to the upstream
	cleaved residue.
	Presently hard-coded for HCV protease
	"""
	ft = FoldTree()
	ft.add_edge(72, 1, -1)
	ft.add_edge(72, 196, -1)
	ft.add_edge(72, 203, 1)
	ft.add_edge(203 ,197, -1)
	ft.add_edge(203 ,207, -1)
	assert ft.check_fold_tree()

	return ft


def make_move_map(pose, peptide_residues, repackable_residues):
	""" 
	Makes a movemap for a protease-peptide system, with all non-peptide 
	residue backbones fixed, and side chains mobile for the peptide and all
	residues in an input selection, which is intended to be the nearby  
	residues, excluding the catalytic ones. 
	"""
	mm = MoveMap()
	for i in peptide_residues:
		mm.set_bb(i, True)
	for i in repackable_residues:
		mm.set_chi(i, True)

	return mm


def make_task_factory(residue_selectors, confine_design=None):
	""" 
	Makes a TaskFactory with operations that leave the mutable residues 
	designable, restricts the nearby residues to repacking, and prevents 
	repacking of other residues.

	Also includes the ability to take in a target residue mutatations text 
	file in the form:
	res_number 		allowed_AAs
	for example:
	138 			KR

	All designable residues not listed in the file are restricted to repacking
	and the listed residues are limited in their design options to those AAs 
	listed.
	"""
	design_set = residue_selectors['mutable']
	repack_set = residue_selectors['packable']
	other_set = residue_selectors['other']

	prevent = PreventRepackingRLT() # No repack, no design
	repack = RestrictToRepackingRLT() # No design

	tf = TaskFactory()
	tf.push_back(OperateOnResidueSubset(prevent, other_set))
	tf.push_back(OperateOnResidueSubset(repack, repack_set))
	# Everything else left designable by default

	# Restricting design further
	if confine_design:
		trms = readfile(confine_design)
		# Converting lines to a dict
		limited_set = \
			{int(line.split()[0]):line.split()[1] for line in trms}

		# Converting residues in the dict to a selector
		res_concatenated = str(list(limited_set.keys())).strip('[]').replace(' ','')
		small_des_set = ResidueIndexSelector(res_concatenated)
		now_only_repack = NotResidueSelector(small_des_set)

		# Making residue selection excluding residues in the file and 
		# restricting them to repacking
		no_longer_designable = AndResidueSelector()
		no_longer_designable.add_residue_selector(design_set)
		no_longer_designable.add_residue_selector(now_only_repack)
		tf.push_back(OperateOnResidueSubset(repack, no_longer_designable))

		# Limiting design on residues in the file
		for res, AAs in list(limited_set.items()):
			designable_res = ResidueIndexSelector(str(res))
			restrict_AAs = RestrictAbsentCanonicalAASRLT()
			restrict_AAs.aas_to_keep(AAs)
			tf.push_back(OperateOnResidueSubset(restrict_AAs, designable_res))

	return tf


def fastrelax(pose, score_function, movemap):
	""" 
	Runs the FastRelax protocol on a pose, using given score function and 
	movemap
	"""
	relax = FastRelax()
	relax.set_scorefxn(score_function)
	relax.set_movemap(movemap)

	relax.apply(pose)
	return pose


def minmover(pose, score_function, movemap):
	""" 
	Runs a gradient-base minimization on a pose, using given score function
	and movemap
	"""	
	min_mover = MinMover()
	min_mover.score_function(score_function)
	min_mover.movemap(movemap)

	min_mover.apply(pose)
	return pose	


def design_pack(pose, score_function, task):
	""" Runs packing mover on a pose, using given score function and task """
	pack = PackRotamersMover(score_function, task)
	pack.apply(pose)
	return pose


def custom_design(pose, score_function, movemap, taskfactory, cycles):
	""" Simple design protocol without fa_rep ramping """

	for i in range(cycles):
		task = taskfactory.create_task_and_apply_taskoperations(pose)
		pose = design_pack(pose, score_function, task)
		pose = minmover(pose, score_function, movemap)

	return pose


def fastdesign(pose, score_function, movemap, taskfactory):
	fd = FastDesign()
	fd.set_scorefxn(score_function)
	fd.set_movemap(movemap)
	fd.set_task_factory(taskfactory)

	fd.apply(pose)
	return pose


######### Evaluation #########################################################
# Weights from the score function REF_2015
ref_wts = \
	{'fa_atr': 1.0, 'fa_rep': 0.55, 'fa_sol': 1.0, 'fa_intra_rep': 0.005, 
	'fa_intra_sol_xover4': 1.0, 'lk_ball_wtd': 1.0, 'fa_elec': 1.0, 
	'pro_close': 1.25, 'hbond_sr_bb': 1.0, 'hbond_lr_bb': 1.0, 
	'hbond_bb_sc': 1.0, 'hbond_sc': 1.0, 'dslf_fa13': 1.25, 
	'atom_pair_constraint': 1.0, 'coordinate_constraint': 1.0, 
	'angle_constraint': 1.0,'dihedral_constraint': 1.0, 'omega': 0.4, 
	'fa_dun': 0.7, 'p_aa_pp': 0.6, 'yhh_planarity': 0.625, 'ref': 1.0, 
	'chainbreak': 1.0, 'rama_prepro': 0.45, 'res_type_constraint': 1.0}


def cleanhead(header):
	""" 
	When using pose.energies(), the lengths of some score terms are long 
	enough that they run into each other. Consequently, the split() method
	doesn't separate them correctly. This function takes a header line and
	cleans up any such clashes identified in ref2015_cst. The reason I bother
	is to accommodate the possibility of cst/no cst, in which case the 
	indexing of the energies table changes.
	"""
	for n, i in enumerate(header):
		# Cleaning up string length clashes
		if i == 'fa_intra_repfa_intra_sol_xo':
			header[n] = 'fa_intra_rep'
			header.insert(n+1, 'fa_intra_sol_xover4')
		if i == 'dslf_fa13atom_pair_constcoordinate_consangle_constraindihedral_constr':
			header[n] = 'dslf_fa13'
			header.insert(n+1, 'atom_pair_constraint')
			header.insert(n+2, 'coordinate_constraint')
			header.insert(n+3, 'angle_constraint')
			header.insert(n+3, 'dihedral_constraint')
		if i == 'rama_preprores_type_constr':
			header[n] = 'rama_prepro'
			header.insert(n+1, 'res_type_constraint')

	return header


def res_scores(pose, residues, score_function):
	""" Gets the residue scores for a given set of residues in a pose """
	score_function(pose)
	pose_energies = str(pose.energies()).split('\n') # Table of residue
		# score components, including a header line, so index matches res

	# Getting scores header
	head_raw = pose_energies[0].split()[1:]
	header = cleanhead(head_raw)

	energies_set = []
	for i in residues:
		res_energies = pose_energies[i].split()[1:]
		res_tot_energy = 0.0
		for j, term in enumerate(header):
			res_tot_energy += (float(res_energies[j]) * ref_wts[term])
		energies_set.append(res_tot_energy)

	set_energy = sum(energies_set)
	return set_energy, energies_set


def move_apart(pose, peptide_start, peptide_end):
	""" Moves the peptide a long distance away from the protease """
	# Making displacement vector
	xyz = xyzVector_double_t()
	xyz.x, xyz.y, xyz.z = [100 for i in range(3)]

	# Moving peptide, atom by atom
	for res in range(peptide_start, peptide_end + 1):
		for atom in range(1, pose.residue(res).natoms() + 1):
			pose.residue(res).set_xyz(atom, pose.residue(res).xyz(atom) + xyz)

	return pose


def score_ddg(pose, tf):
	"""
	Gets a score for the ddG of peptide binding to the protease. This is 
	achieved by taking the peptide and moving each atom a set arbitrary length 
	that is large enough to be far from the protease, then repacking the side
	chains of both the peptide and the mutable protease residues. This 
	function does not take a scorefunction as an input, scoring instead with 
	the default function to ignore the effects of constraints.
	"""
	# Making a new pose to avoid messing up the input
	ddg_pose = Pose()
	ddg_pose.assign(pose)

	sf = get_fa_scorefxn()

	# Score when docked
	dock_score = sf(ddg_pose)

	# Score when separate
	pt = tf.create_task_and_apply_taskoperations(pose)
	pt.restrict_to_repacking()
	ddg_pose = move_apart(ddg_pose, 197, 207)
	ddg_pose = design_pack(ddg_pose, sf, pt)
	split_score = sf(ddg_pose)

	ddg = dock_score - split_score
	return round(ddg, 3)


def ident_mutations(start_pose, end_pose, residues, start_set, end_set):
	"""
	Compares the sequences of a starting pose and ending pose at specified 
	residues and identifies differences. Returns a string listing changes in 
	the format of ANB, where A is the starting residue, N is the residue 
	number, and B is the ending residue.
	"""
	template = '{:6s}{:6s}{:6s}{:12s}'
	mutations = ''
	mutations_present = False
	for i in residues:
		res_mut_info = [i]
		start_res = start_pose.residue(i).name1()
		res_mut_info.append(start_res)
		end_res = end_pose.residue(i).name1()
		if start_res != end_res:
			mutations_present = True
			res_mut_info.append(end_res)

			r_i = residues.index(i)
			e_dif = round(end_set[r_i] - start_set[r_i], 3)
			res_mut_info.append(e_dif)

		else:
			res_mut_info += ['NO CHANGE', '']

		mutations += template.format(*[str(i) for i in res_mut_info])

	if mutations_present:
		return mutations.lstrip(',')
	else:
		return "NONE"


def read_pose_scores(pose):
	""" Gets weighted scores for a full pose, returns as list of strings """
	pose_scores = []
	energies = pose.energies().total_energies()
	energies_as_list = [i.strip('( )') for i in str(energies).split(') (')]

	for e in energies_as_list:
		term, unweighted_val = e.split()
		term = term.replace(';', '')
		if term in ref_wts:
			weighted_val = float(unweighted_val) * ref_wts[term]
			pose_scores.append(': '.join([term,str(weighted_val)]))

	return pose_scores


######### Design and evaluation protocol #####################################
class mutation_collection:
	"""
	This class object will include information on the mutations in decoys 
	produced by design_protease.py

	Decided against collecting all data this way because it doesn't work well 
	with parallel processes on larger sets.
	"""
	def __init__(self, threaded_pose, selector_dict):
		self.threaded_pose = threaded_pose
		self.selectors = selector_dict
		self.selectors_to_list()

		self.decor_identifiers = []
		self.relaxed_decoys = []
		self.designed_decoys = []

	def selectors_to_list(self):
		""" Converts a selector output vector to lists of selected residues """
		for set_name in self.selectors:
			selector = self.selectors[set_name]
			selection_list = selector_to_list(self.threaded_pose, selector)
			
			setattr(self, set_name + '_residues', selection_list)

	def read_sequence_name(self, name):
		""" 
		Reads name of file generated earlier in this script by quick_thread 
		and determines whether it is a cleaved or uncleaved sequence and what
		the sequence is.
		"""
		breakup = name.split('_')
		print(breakup)

		# Determining whether sequence is cleaved
		cleavage = breakup[0]
		assert cleavage in ['cleaved', 'uncleaved']
		if cleavage == 'cleaved':
			self.cleaved = True
		else:
			self.cleaved = False

		# Determining sequence
		self.sequence = breakup[-2]
		self.recognition_sequence = self.sequence[1:7]


def set_design(pdb, residue_selectors, args):
	"""
	Uses the job distributor to output a set of proteases designed for 
	compatability with a threaded peptide sequence. Outputs a provided number
	of decoys into a directory that is also a required input. Will relax the 
	pdb, then run 20 rounds of design/repacking plus minimization. For all 
	movers, only the residues in the peptide and those within the input list
	are repackable, and only those in the input list are designable. For the 
	relax, the peptide backbone is flexible, and constraints are applied.
	"""
	pose = apply_constraints(pose_from_pdb(pdb))
	mc = mutation_collection(pose, residue_selectors)

	sf = create_score_function('ref2015_cst')
	ft = make_fold_tree() # Hard-coded for HCV protease
	pose.fold_tree(ft) # Improve sampling efficiency
	mm = make_move_map(pose, mc.peptide_residues, mc.movemap_pack_residues)
	tf = make_task_factory(residue_selectors, args.target_res_mutation_file)

	dec_name = pdb.replace('.pdb.gz', '_designed')
	jd = PyJobDistributor(dec_name, args.number_decoys, sf)

	while not jd.job_complete:
		pp = Pose()
		pp.assign(pose)

		# Relaxing
		print('relaxing')
		relax_name = jd.current_name.replace('designed', 'relaxed')
		pp = fastrelax(pp, sf, mm)
		relaxed_struc = Pose()
		relaxed_struc.assign(pp)
		relaxed_struc.dump_pdb(relax_name)

		# Doing design, default is FastDesign
		print('designing')
		if args.method == 'custom':
			pp = custom_design(pp, sf, mm, tf, 20)
		else:
			pp = fastdesign(pp, sf, mm, tf)		

		# Getting residue scores, ddG, and mutations list
		score_change = sf(pp) - sf(relaxed_struc)
		relax_res_E_sum, relax_res_E_set = \
			res_scores(relaxed_struc, mc.mutable_residues, sf)
		des_res_E_sum, des_res_E_set = res_scores(pp, mc.mutable_residues, sf)
		res_score_change = des_res_E_sum - relax_res_E_sum
		ddg = score_ddg(pp, tf)
		ddg_change = ddg - score_ddg(relaxed_struc, tf)
		pro_mut = ident_mutations(pose, pp, mc.mutable_residues, 
									relax_res_E_set, des_res_E_set)

		# Making line to add to fasc file
		scores = [score_change, des_res_E_sum, res_score_change, ddg, ddg_change, pro_mut]
		temp = "score_change {}   residue_scores: {}   residue_score_change: {}   ddG: {}   ddG_change: {}   mutations: {}"
		score_text = temp.format(*[str(i) for i in scores])
		print(score_text, '\n')
		jd.additional_decoy_info = score_text

		jd.output_decoy(pp)

######### Main ###############################################################
def main():
	# Getting user inputs
	args = parse_args()

	# Initializing PyRosetta
	ros_opts = init_opts(cst_file=args.constraints)
	init(options=ros_opts)

	# Destination folder for PDB files
	dir_nam = args.out_dir
	if not isdir(dir_nam):
		makedirs(dir_nam)

	# Reading inputs for peptide sequences
	cut_seq = get_seq_list(args.cut_peptide_sequence)
	uncut_seq = get_seq_list(args.uncut_peptide_sequence)

	# Creating threaded structures
	pose = pose_from_pdb(args.start_struct)
	make = args.thread
	if make:
		print("Threading")
	t_structs = quick_thread(dir_nam, pose, cut_seq, cleaved=True, make=make)
	t_structs += quick_thread(dir_nam, pose, uncut_seq, make=make)

	# Making residue selectors
	residue_selectors = select_residues(design_peptide=args.design_peptide)

	# Doing design on threaded models
	for struc in t_structs:
		set_design(struc, residue_selectors, args)


if __name__ == '__main__':
	main()