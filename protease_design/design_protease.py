#!/usr/bin/python
""" 
PyRosetta4, Python 3.5
Joseph Lubin, 2019

Pipeline for rapidly modeling protease-substrate combinations, and using 
FastRelax and FastDesign to explore potentially better interacting variants

It is assumed in the program that the input PDB structure will have a set of
enzdes constraint comments at the beginning of the document, and that the 
protease is chain A and the substrate is chain B.

Sample command:
python design_protease.py -s HCV.pdb -od F2R2_des_pep -name LY104_F2R2_test_pep  
-seq DVDAR -site 198 -ps "198-202" -cons ly104.cst -nd 100 -mm 138 I -mm 170 Q  
-mm 171 S -mm 173 I -mm 175 K -mm 183 R -pep_only
"""
from __future__ import print_function # For compatability with Python 2.7
import argparse
from os import makedirs
from os.path import basename, isdir, isfile, join
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
	PreventRepackingRLT, RestrictAbsentCanonicalAASRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.select.residue_selector import \
	AndResidueSelector, ChainSelector, InterGroupInterfaceByVectorSelector,\
	NotResidueSelector, OrResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.teaching import SimpleThreadingMover
from random import randint
from sys import exit

def parse_args():
	info = "Design a protease around a peptide sequence"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-s", "--start_struct", required=True,
		default='start_proteases/HCV.pdb', help="Pick starting PDB")
	parser.add_argument("-od", "--out_dir", required=True,
		help="Name an output directory for decoys")
	parser.add_argument("-name", "--name", type=str,
		help="How would you like to name your outputs? \
		(Default will use the name of the input PDB file.)")
	parser.add_argument("-seq", "--sequence", required=True, type=str,
		help="What substrate sequence do you want to thread")
	parser.add_argument("-site", "--subst_site", required=True, type=int,
		help="Where in the pose should the substrate sequence begin \
		threading? (A 5-letter substitution ending with 201 should start at \
		197).")
	parser.add_argument("-cr", "--cat_res", type=int, nargs='+', 
		default=[72, 96, 154], help="The catalytic residues of the protease, \
		excluded from design (defaults are 72, 96, and 154, for HCV)")
	parser.add_argument("-ps", "--pep_subset", type=str, default=None, 
		help='Select the subset of the peptide around which to design, as a \
		string of "first_res-last_res". (Ex: "198-202") Otherwise, design \
		will be performed around the full peptide. Both options that include \
		peptide design will only design the subset.')
	parser.add_argument("-cons", "--constraints", type=str, 
		default='ly104.cst', help="Pick constraints file")
	parser.add_argument("-des_pep", "--design_peptide", action="store_true", 
		help="Option to allow design on the regognition region of the peptide.")
	parser.add_argument("-pep_only", "--design_only_peptide", 
		action="store_true", help="Option to allow design only of the \
		peptide, excluding the surrounding protease.")
	parser.add_argument("-no_design", "--no_design", 
		action="store_true", help="Option to just modify the peptide and \
		protease to a desired sequence and relax.")
	parser.add_argument("-hbn", "--use_hb_net", action="store_true", 
		help="Option to include HBnet score term in design.")
	parser.add_argument("-n", "--number_decoys", type=int, default=10, 
		help="How many decoys should be made? (Default is 10.)")
	parser.add_argument("-mm", "--mutations", nargs=2, action='append',
		help="Manually input mutations in the format [site] [one-letter res]. \
		Accepts multiple uses. (Ex: -mm 138 I -mm 183 R) Note, if you intend \
		to change the catalytic residues, you must edit the PDB's enzdes \
		comments as well, or applying constraints won't work properly.")
	parser.add_argument("-trmf", "--target_res_mutation_file", type=str, 
		default=None, help="File of format: 'res_number    design_options' \
		to specify further design restrictions. Useful when picking between \
		identified favorable mutations.")
	parser.add_argument("-test", "--test_mode", action="store_true", 
		help="For debugging: test protocol, exiting before generating decoys.")
	args = parser.parse_args()
	return args


def init_opts(cst_file='ly104.cst'):
	""" Produces a list of init options for PyRosetta, including cst file """
	ros_opts = '-ex1 -ex2  -use_input_sc -flip_HNQ'
	ros_opts += ' -mute all -enzdes::cstfile {}'.format(cst_file)
	ros_opts += ' -cst_fa_weight 1.0 -run:preserve_header -out:pdb_gz'
	return ros_opts


def readfile(file_name):
	""" Opens a file in read-mode and returns a list of the text lines """
	with open(file_name, 'r') as r:
		lines = r.readlines()

	return lines

######### Threading ##########################################################

def input_manual_mutations(pose, mutations):
	"""
	Allows the manual input of specific mutations into a starting pose
	Takes input pose and a list of mutations, with each mutation as a two-item
	list, [site, single-letter residue name].
	"""
	mutated_pose = Pose(pose)
	print('\nInputting {} mutations:'.format(str(len(mutations))))
	for m in mutations:
		print('\t{}{}{}'.format(pose.residue(int(m[0])).name1(), m[0], m[1]))
		mmover = SimpleThreadingMover(m[1].upper(), int(m[0]))
		mmover.apply(mutated_pose)

	return mutated_pose


def random_aa(length):
	""" 
	Returns a string of random 1-letter amino acid names from the cannonical 
	20, to a specified length.
	"""

	aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
				'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	
	aa_string = ""

	for aa in range(length):
		rand_index = randint(0,19)
		aa_string += aa_list[rand_index]

	return aa_string


def thread_substrate(destination, name, pose, sequence, peptide_start):
	"""
	Creates a threaded PDB from a given destination, name, and sequence and 
	returns a threaded pose. Requires peptide start site to know where to begin
	threading.
	"""
	print('\nThreading {} at site {}'.format(sequence, str(peptide_start)))

	# Making threaded pose
	threaded_pose = Pose(pose)
	tm = SimpleThreadingMover(sequence.upper(), peptide_start)
	tm.apply(threaded_pose)

	# Outputting threaded pose
	threaded_out = '_'.join([name, 'threaded.pdb'])
	pdb_name = join(destination, threaded_out)
	threaded_pose.dump_pdb(pdb_name)

	print('Saved as {}'.format(pdb_name))
	return threaded_pose

######### Residue selection ##################################################

def selector_intersection(*selectors):
	""" Returns the intersection of any set of selectors """
	intersect_selection = AndResidueSelector()
	for s in selectors:
		intersect_selection.add_residue_selector(s)

	return intersect_selection


def selector_union(*selectors):
	""" Returns the intersection of any set of selectors """
	union_selection = OrResidueSelector()
	for s in selectors:
		union_selection.add_residue_selector(s)

	return union_selection


def mutable_residues_selector(protease_selection, peptide_selection,
	catalytic_selection, design_peptide=False):
	"""
	Selects the residues in a shell around the peptide using the 
	InterGroupInterfaceByVectorSelector residue selector
	"""
	# Making protease shell selector (includes peptide)
	first_shell_select = InterGroupInterfaceByVectorSelector()
	first_shell_select.group1_selector(protease_selection) 
	first_shell_select.group2_selector(peptide_selection) 
	first_shell_select.nearby_atom_cut(8)
	first_shell_select.vector_dist_cut(10)

	# Excluding the catalytic residues, peptide (if not designed)
	not_cats_sel = NotResidueSelector(catalytic_selection)
	if design_peptide:
		mutable_selection = selector_intersection(
			not_cats_sel, first_shell_select)
	else:
		mutable_selection = selector_intersection(
			not_cats_sel, first_shell_select, protease_selection)

	return mutable_selection


def packable_residues_selector(
	peptide_selection, mutable_selection, catalytic_selection):
	"""
	Selects the shell of neighbor residues to repack. Packable set should not
	include the mutable set, since the action is RestrictToRepacking.
	"""
	# Making negative selections for mutable and catalytic
	not_mutable = NotResidueSelector(mutable_selection)
	not_catalytic = NotResidueSelector(catalytic_selection)

	# Selecting residues near mutable shell
	near_mutable = InterGroupInterfaceByVectorSelector()
	near_mutable.group1_selector(not_mutable)
	near_mutable.group2_selector(mutable_selection)
	near_mutable.nearby_atom_cut(8)
	near_mutable.vector_dist_cut(10)

	# Selecting residues near the peptide, with wider range for BB mobility
	near_pep = InterGroupInterfaceByVectorSelector()
	near_pep.group1_selector(not_mutable)
	near_pep.group2_selector(peptide_selection)
	near_pep.nearby_atom_cut(10)
	near_pep.vector_dist_cut(12)

	# Combining selections for peptide and near peptide and near mutable
	inclusive_packable = selector_union(
		near_mutable, near_pep, peptide_selection, ChainSelector('B')) ##################Chain B hacky
	
	# Setting up exclusion of catalytic and mutable residues
	exclusive_packable = selector_intersection(
		inclusive_packable, not_mutable, not_catalytic)

	return exclusive_packable


def select_residues(cat_res, peptide_subset, design_peptide=False, 
	design_protease=True):
	""" 
	Makes residue selectors for protease sections. Requires manual input for 
	which residues are catalytic and whether only part of the peptide should be 
	selected. Options for whether the peptide is designable (false by default) 
	and whether the protease is designable (true by default). Assumes that the 
	protease is chain A and the peptide is chain B.
	"""
	residue_selectors = {}

	# Catalytic residues. ResidueIndexSelector needs a string, not a list.
	cats_as_str = ','.join([str(i) for i in cat_res]) 
	catalytic = ResidueIndexSelector(cats_as_str)
	residue_selectors['catalytic'] = catalytic

	# Protease residues. Protease assumed to be chain A
	protease = ChainSelector("A")
	residue_selectors['protease'] = protease
	
	# Peptide residues. Peptide assumed to be chain B, unless range specified
	if peptide_subset:
		peptide = ResidueIndexSelector(peptide_subset)
	else:
		peptide = ChainSelector("B")
	residue_selectors['peptide'] = peptide

	# Designable residues. May include protease and peptide, just one, or none
	if design_protease:
		mutable = mutable_residues_selector(protease, peptide,
			catalytic, design_peptide)
	elif design_peptide:
		mutable = peptide
	else: # Neither protease not peptide designable
		mutable = selector_intersection(protease, peptide) # Empty set
	residue_selectors['mutable'] = mutable

	# Packable residues. Centered around the peptide and designable set
	packable = packable_residues_selector(peptide, mutable, catalytic)
	residue_selectors['packable'] = packable
	
	# Immobile residues. Catalytic residues and everything that isn't mutable 
	# or packable
	immobile = NotResidueSelector(selector_union(mutable, packable))
	residue_selectors['immobile'] = immobile

	return residue_selectors


def selector_to_list(pose, selector):
	""" Converts a selector output vector to a list of selected residues """
	selection_vector = selector.apply(pose)
	selection_list = []
	for i in range(len(selection_vector)): 
		if selection_vector[i+1]==1:
			selection_list.append(i+1)

	return selection_list 

######### Setup ##############################################################

def apply_constraints(pose):
	""" Applies the constraints form the input CST file to a pose """
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)
	return pose


def make_move_map(pose, selectors):
	""" 
	Makes a movemap for a protease-peptide system, with all non-peptide 
	residue backbones fixed, and side chains mobile for all packable and 
	designable residues. 

	Takes a dict of selectors generated by select_residues.
	"""
	mm = MoveMap()

	# Mobile backbone for peptide
	for i in selector_to_list(pose, selectors['peptide']):
		mm.set_bb(i, True)
	
	# Mobile side chains for mutable and packable residues
	for i in selector_to_list(pose, selectors['mutable']):
		mm.set_chi(i, True)
	for i in selector_to_list(pose, selectors['packable']):
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
	mutable_set = residue_selectors['mutable']
	repack_set = residue_selectors['packable']
	immobile_set = residue_selectors['immobile']

	prevent = PreventRepackingRLT() # No repack, no design
	repack = RestrictToRepackingRLT() # No design

	tf = TaskFactory()
	tf.push_back(IncludeCurrent())
	tf.push_back(ExtraRotamers(0, 1, 1))
	tf.push_back(ExtraRotamers(0, 2, 1))
	tf.push_back(OperateOnResidueSubset(prevent, immobile_set))
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
		no_longer_designable.add_residue_selector(mutable_set)
		no_longer_designable.add_residue_selector(now_only_repack)
		tf.push_back(OperateOnResidueSubset(repack, no_longer_designable))

		# Limiting design on residues in the file
		for res, AAs in list(limited_set.items()):
			designable_res = ResidueIndexSelector(str(res))
			restrict_AAs = RestrictAbsentCanonicalAASRLT()
			restrict_AAs.aas_to_keep(AAs)
			tf.push_back(OperateOnResidueSubset(restrict_AAs, designable_res))

	return tf


def get_score_function(constraints=True, hbnet=False):
	""" Returns either default or weighted REF2015 with or without hbnet """
	# Picking between constraints and not
	if constraints:
		sf = create_score_function('ref2015_cst')
	else:
		sf = create_score_function('ref2015')

	# Optionally adding in hbnet
	if hbnet:
		sf.set_weight(ScoreType.hbnet, 1)

	return sf

######### Design Protocols ###################################################

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


def fastdesign(pose, score_function, movemap, taskfactory):
	fd = FastDesign()
	fd.set_scorefxn(score_function)
	fd.set_movemap(movemap)
	fd.set_task_factory(taskfactory)

	fd.apply(pose)
	return pose


def jd_design(name, decoy_count, pose, score_function, movemap, task_factory, 
	do_design=True):
	""" Runs job distributor with relax and design protocols """
	print('\n')

	jd = PyJobDistributor(name, decoy_count, score_function)
	while not jd.job_complete:
		pp = Pose(pose)

		# Relaxing
		print('Relaxing...')
		pp = fastrelax(pp, score_function, movemap)

		if do_design:
			relax_name = jd.current_name.replace('designed', 'relaxed')
			pp.dump_pdb(relax_name)

		# Doing design and outputting decoy
		print('Designing...')
		pp = fastdesign(pp, score_function, movemap, task_factory)	

		print('Complete\n')
		jd.output_decoy(pp)

	return

######### Main ###############################################################

def test_and_exit(args, residue_selectors, pose, name):
	""" Prints info then exits """
	print('\n\nArgs:')
	print(args)
	print('\nSelectors:')
	for k, v in residue_selectors.items():
		print('\t',k)
		print('\t',selector_to_list(pose,v))
	print('\nSequence')
	print(pose.sequence())
	print('\nName:')
	print(name)
	print('\n')

	exit()


def main(args):
	# Initializing PyRosetta
	ros_opts = init_opts(cst_file=args.constraints)
	init(options=ros_opts)

	# Destination folder for PDB files
	dir_name = args.out_dir
	if not isdir(dir_name):
		print('\nMaking directory: {}'.format(dir_name))
		makedirs(dir_name)

	# Getting name for outputs
	if args.name:
		out_name = args.name
	else:
		out_name = basename(args.start_struct)
		# strip out .pdb or .pdb.gz extension
		out_name = out_name.replace('.pdb', '').replace('.gz', '')	

	# Preparing pose, with manual mutations, substrate threading, constraints
	pose = pose_from_pdb(args.start_struct)

	if args.mutations: # Making manually input mutations
		pose = input_manual_mutations(pose, args.mutations)
	
	# If the peptide is to be designed, starting with a random peptide sequence
	if args.design_peptide or args.design_only_peptide: ############################# put new substrate in for every run?
		seq_to_thread = random_aa(len(args.sequence))
		design_peptide = True
	else:
		seq_to_thread = args.sequence
		design_peptide = False

	pose = thread_substrate(
		dir_name, out_name, pose, seq_to_thread, args.subst_site)

	pose = apply_constraints(pose)

	# Making residue selectors
	design=True
	if args.no_design:
		design = False

	design_protease = True
	if args.design_only_peptide or args.no_design:
		design_protease=False

	residue_selectors = select_residues(args.cat_res, args.pep_subset, 
			design_peptide=design_peptide, design_protease=design_protease)

	# Creating score function, movemap, and taskfactory for design
	sf = get_score_function(constraints=True, hbnet=args.use_hb_net)
	mm = make_move_map(pose, residue_selectors)
	tf = make_task_factory(residue_selectors, args.target_res_mutation_file)

	# Running relax and design protocol
	dec_name = join(dir_name, out_name)
	if design:
		dec_name += '_designed'

	if args.test_mode:
		test_and_exit(args, residue_selectors, pose, dec_name)

	jd_design(dec_name, args.number_decoys, pose, sf, mm, tf, do_design=design)


if __name__ == '__main__':
	args = parse_args()
	main(args)