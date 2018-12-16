#!/usr/bin/python
""" 
PyRosetta4, Python 2.7
Joseph Lubin, 2017
"""
import argparse
from math import sqrt
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta.teaching import SimpleThreadingMover
from pyrosetta.teaching import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import parse_resfile
from pyrosetta.rosetta.protocols.enzdes import AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW
from glob import glob
from os import makedirs
from os.path import basename, isdir, join

# 72 and 154
def parse_args():
	info = "Design a protease around a peptide sequence"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-s", "--start_struct", required=True,
		help="Pick starting PDB")
	parser.add_argument("-cseq", "--cut_peptide_sequence", type=str, 
		action='append', help="List cleaved peptide sequences or provide \
		list file")
	parser.add_argument("-useq", "--uncut_peptide_sequence", type=str, 
		action='append', help="List uncleaved peptide sequences or provide \
		list file")
	parser.add_argument("-cr", "--cat_res", type=str, nargs='+', 
		default=[72, 154], help="The catalytic residues of the protease, \
		excluded from design (defaults are 72 and 154, for HCV)")
	parser.add_argument("-cst", "--constraints", type=str,
		help="Pick constraints file")
	parser.add_argument("-rad", "--radius_design", type=int, default=8, 
		help="Cutoff for designable residues (Angstroms from peptide, \
		default is 8A)")
	parser.add_argument("-th", "--thread", action = "store_true", 
		help="Option to create threaded models. Use the first time running.")
	#parser.add_argument("-re", "--rerun", type=str, 
	#	help="If the protocol has already been run, in what directory \
	#	are the _best.pdb structures stored?")
	parser.add_argument("-res", "--resfile", type=str,
		help="Pick resfile for design")
	args = parser.parse_args()
	return args


def init_opts(cst_file):
	""" Produces a list of init options for PyRosetta, including cst file """
	ros_opts = '-mute core -mute protocols -mute basic'
	ros_opts += ' -enzdes::cstfile ' + cst_file
	ros_opts += ' -cst_fa_weight 1.0 -run:preserve_header -out:pdb_gz'
	return ros_opts


def res_ca_cords(pose, res):
	""" Returns the x,y,z coordinates of the A-carbon of a given residue"""
	res_coords = []
	for i in range(3):
		res_coords.append(pose.residue(res).xyz('CA')[i])

	return res_coords


def point_dist(point_1, point_2):
	""" Given two sets of coordinates, determines distance between them """
	sum_difs_squared = 0
	for i in range(3):
		sum_difs_squared += (point_2[i] - point_1[i]) ** 2

	return sqrt(sum_difs_squared)


def res_to_design(pdb, radius, cat_res):
	""" 
	Determines the chain in the PDB that is the peptide, assuming that the 
	peptide is smallest. Determines the coordinates of all a-carbons in the 
	system, then checks the distances to each atom in the peptide. If the 
	atom is within the cutoff radius, it is added to the list of designable 
	residues. Returns a list of mutable residues.
	"""
	pose = pose_from_pdb(pdb)
	chains = pose.split_by_chain()

	# Determining peptide chain
	pep_chain_no = 1
	for i in range(2, len(chains) + 1):
		if len(chains[i]) < len(chains[pep_chain_no]):
			pep_chain_no = i
	pep_chain = chains[pep_chain_no]
	pep_len = pep_chain.total_residue()
	chains.pop(pep_chain_no) # removes peptide chain from chain list

	# Getting peptide residue CA coordinates:
	pep_coords = []
	for i in range(1, pep_len + 1):
		pep_coords.append(res_ca_cords(pep_chain, i))

	# Populating the list of designable residues
	mutable_residues = []
	for chain in chains:
		# Note that chains now excludes the peptide chain
		for res in range(1, pep_len + 1):
			if res in cat_res:
				# Exclude the catalytic residues from the designable list
				continue
			a_cords = res_ca_cords(pose, res)
			for i in pep_coords:
				# If any res is within radius of any pep res, add to the list
				if point_dist(a_cords, i) <= radius:
					mutable_residues.append(res)
					break

	return mutable_residues 


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
			with open(inp, 'r') as t:
				lis = t.readlines()
			if len(lis) == 1:
				# If all sequences are listed horizontally on one line
				# rather than one per line, rearrange
				lis = lis[0].split()

			for i in lis:
				pep_sequences.append(i.strip())

		else:
			# Sequence was typed directly into the argument
			pep_sequences.append(inp.strip())

	return pep_sequences
		

def thread_seq(pose, pep_start, pep_length, seq):
	""" Thread a new sequence in for the peptide. """
	tm = SimpleThreadingMover(seq, pep_start)
	tm.apply(pose)
	return pose


def quick_thread(destination, pdb, sequences, cleaved=False, make=False):
	""" 
	Threads a set of sequences onto the peptide portion of the given PDB file,
	outputting a threaded PDB file for each sequence.
	Function is presently hard-coded for this application.
	"""
	pose = pose_from_pdb(pdb)

	thread_files = []

	for seq in sequences:
		# Naming model
		if cleaved:
			pdbname = 'cleaved_ly104_wt_' + seq + '.pdb.gz'
		else:
			pdbname = 'uncleaved_ly104_wt_' + seq + '.pdb.gz'
		out_name = join(destination, pdbname)
		thread_files.append(out_name)

		if make:
			# Threading peptide sequences
			threaded_pose = Pose()
			threaded_pose.assign(pose)
			threaded_pose = thread_seq(threaded_pose, 197, 11, seq)
			threaded_pose.dump_pdb(out_name)

	return thread_files


def apply_constraints(pose):
	""" Applies the constraints form the input CST file to a pose """
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)
	return pose


def make_move_map(pep_start, pep_end, near_res):
	""" 
	Makes a movemap for a protease-peptide system, with all non-peptide 
	residue backbones fixed, and side chains mobile for the peptide and all
	residues in an input list, which is intended to be the nearby residues 
	(8A by default), excluding the catalytic ones. 
	"""
	mm = MoveMap()
	mm.set_bb_true_range(pep_start,pep_end)
	mm.set_chi_true_range(pep_start,pep_end)
	for i in near_res:
		mm.set_chi(i, True)
	
	return mm


def make_pack_task(pose, near_res=[], resfile=None, design=True):
	"""
	Makes a packer task for a given pose. Has different behaviors for peptide,
	protease, or peptide + protease poses (of type determined by the number of 
	residues). For a peptide (<15 residues), restricts to repacking. For a 
	protease (<200 residues), restricts to repacking and fixes backbone for all
	residues except those in a provided list. For a full system, reads an 
	input resfile.
	"""
	size = pose.total_residue()

	task = standard_packer_task(pose)
	if size < 15:
		# Packer for peptide
		task.restrict_to_repacking()

	elif size < 200:
		# Packer for protease
		task.restrict_to_repacking()
		task.temporarily_fix_everything()
		for i in near_res:
			task.temporarily_set_pack_residue(i, True)

	else:
		# Packer for protease + peptide\
		parse_resfile(pose, task, resfile)
		if design == False:
			task.restrict_to_repacking()


	return task


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


def res_scores(pose, residues, score_function):
	""" Gets the residue scores for a given set of residues in a pose """
	sf = score_function
	sf(pose)
	pose_energies = str(pose.energies()).split('\n') # Table of residue
		# score components, including a header line, so index matches res
	set_energy = 0
	for i in residues:
		res_energies = pose_energies[i].split()[1:]
		res_tot_energy = sum([float(j) for j in res_energies])
		set_energy += res_tot_energy

	return set_energy


def score_ddg(pose, near_res):
	"""
	Gets a score for the ddG of peptide binding to the protease. This is 
	achieved by taking the peptide and moving each atom a set arbitrary length 
	that is large enough to be far from the protease, then repacking the side
	chains of both the peptide and the mutable protease residues. This 
	function does not take a scorefunction as an input, scoring instead with 
	the default function to ignore the effects of constraints.
	"""


def set_design(pdb, score_function, des_res, num_decoys, resfile):
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
	sf = score_function	
	mm = make_move_map(197, 207, des_res) # for relax and minimization
	pt = make_pack_task(pose, near_res=des_res, resfile=resfile) # for packer
	print mm
	print pt

	dec_name = pdb.replace('.pdb.gz', '_relaxed')
	jd = PyJobDistributor(dec_name, num_decoys, sf)

	while not jd.job_complete:
		# Relaxing
		pp = Pose()
		pp.assign(pose)
		pp = fastrelax(pp, sf, mm)

		# Doing design
		for i in range(20):
			pp = design_pack(pp, sf, pt)
			pp = minmover(pp, sf, mm)

		# Getting residue scores
		prot_res_e = res_scores(pose, des_res, sf)
		pep_res_e = res_scores(pose, range(197, 208), sf)



		jd.output_decoy(pp)


def ddg_score(pose, near_res, score_function):
	"""
	Splits the pose into the peptide and the protease, then repacks the 
	rotamers that were flexible in the initial relax/design step and
	rescores to calculate the ddG.
	"""
	# Splitting pose
	chains = pose.split_by_chain()
	pro_pose = Pose()
	pro_pose.assign(chains[1])
	pep_pose = Pose()
	pep_pose.assign(chains[2])

	sf = score_function

	# Packer tasks
	pro_task = make_pack_task(pro_pose, near_res=near_res)
	pep_task = make_pack_task(pep_pose)

	# Applying movers
	pro_pack = PackRotamersMover(sf, pro_task)
	pro_pack.apply(pro_pose)
	pep_pack = PackRotamersMover(sf, pep_task)
	pep_pack.apply(pep_pose)

	# Scoring
	dock_score = sf(pose)
	pro_score = sf(pro_pose)
	pep_score = sf(pep_pose)
	split_score = pro_score + pep_score 
	ddg = dock_score - split_score

	return [dock_score, split_score, ddg] 


def write_report(pose, seq, des_res, sf, head=False):
	"""	Creates a line for a report file with residue and ddG scores """
	template = '{:15s}' * 6 + '\n'
	head_line = ['sequence', 'prot_res_E', 'pep_res_E', 'docked_score', 
				'split_score', 'ddG']

	if head == True:
		return template.format(*head_line)

	else:
		# Getting scores
		prot_res_e = res_scores(pose, des_res, sf)
		pep_res_e = res_scores(pose, range(197, 208), sf)
		dock_split_ddg = ddg_score(pose, des_res, sf)

		# Generating line
		line = [seq]
		line.append(str(prot_res_e))
		line.append(str(pep_res_e))
		line += [str(i) for i in dock_split_ddg]
		return template.format(*line)


def main():
	# Getting user inputs
	args = parse_args()

	# Initializing PyRosetta
	ros_opts = init_opts(args.constraints)
	init(options=ros_opts)

	# Score function
	sf = create_score_function('ref2015_cst')

	# Destination folder for PDB files
	pdb = args.start_struct
	dir_nam = 'ly104_design_decoys'
	if not isdir(dir_nam):
		makedirs(dir_nam)

	# Reading inputs for peptide sequences
	cut_seq = get_seq_list(args.cut_peptide_sequence)
	uncut_seq = get_seq_list(args.uncut_peptide_sequence)

	# Determining peptide part of PDB file, residues near peptide
	des_res = res_to_design(pdb, args.radius_design, args.cat_res)

	# Creating threaded structures
	make = False
	if args.thread:
		make = True
	t_structs = quick_thread(dir_nam, pdb, cut_seq, cleaved=True, make=make)
	t_structs += quick_thread(dir_nam, pdb, uncut_seq, make=make)
	print t_structs

	"""
	# Doing design on threaded models
	for struc in t_structs:


	thread_and_relax(pdb, pep_sequences, dir_nam, des_res, sf)

	# Getting list of best structures for ddG scoring
	best_structures = get_best_structs(dir_nam)

	# Analyzing best models
	report_name = pdb.replace('.pdb', '_relax_report.txt')
	with open(report_name, 'w') as r:
		for i in range(len(pep_sequences)):
			if i == 0:
				# Writing header
				line = write_report('', '', '', '', head=True)
				r.write(line)
			# Writing report line
			line = write_report(best_structures[i], pep_sequences[i],
								des_res, sf, head=False)
			r.write(line)
	"""


if __name__ == '__main__':
	main()