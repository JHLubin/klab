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
from pyrosetta.rosetta.protocols.enzdes import AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW
from glob import glob
from os import makedirs
from os.path import basename, isdir, join

res_names = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 
			'G':'GLY','H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 
			'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', 
			'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}

block = "REMARK 666 MATCH TEMPLATE A ASP 96 MATCH MOTIF A HIS 72  1  1\n" +\
		"REMARK 666 MATCH TEMPLATE A HIS 72 MATCH MOTIF A SER 154  2  1\n" +\
		"REMARK 666 MATCH TEMPLATE A SER 154 MATCH MOTIF B CYS 7  3  1\n" +\
		"REMARK 666 MATCH TEMPLATE A HIS 72 MATCH MOTIF B SER 8  4  1\n"

def parse_args():
	info = "Design a protease around a peptide sequence"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-s", "--start_struct", required=True,
		help="Pick starting PDB")
	parser.add_argument("-seq", "--peptide_sequence", type=str, 
		action='append', help="List peptide sequences or provide list file")
	parser.add_argument("-cr", "--cat_res", type=str, nargs='+', 
		default=[72, 154], help="The catalytic residues of the protease, \
		excluded from design (defaults are 72 and 154, for HCV)")
	parser.add_argument("-cst", "--constraints", type=str,
		help="Pick constraints file")
	parser.add_argument("-rad", "--radius_design", type=int, default=8, 
		help="Cutoff for designable residues (Angstroms from peptide, \
		default is 8A)")
	parser.add_argument("-re", "--rerun", type=str, 
		help="If the protocol has already been run, in what directory \
		are the _best.pdb structures stored?")
	#parser.add_argument("-res", "--resfile", type=str,
	#	help="Pick resfile, if available")
	args = parser.parse_args()
	return args


def init_opts(cst_file):
	""" Produces a list of init options for PyRosetta, including cst file """
	ros_opts = '-mute core -mute protocols -mute basic -enzdes::cstfile '
	ros_opts += cst_file
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
	residues. Returns the peptide chain starting residue and length, and a 
	list of mutable residues.
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

	# Getting residue number of peptide start
	pep_start = 1
	for i in range(1,pep_chain_no):
		pep_start += chains[i].total_residue()

	# Getting peptide residue CA coordinates:
	pep_coords = []
	for i in range(1, pep_chain.total_residue() + 1):
		pep_coords.append(res_ca_cords(pep_chain, i))

	# Populating the list of designable residues
	mutable_residues = []
	for chain in chains:
		# Note that chains now excludes the peptide chain
		for res in range(1, chain.total_residue() + 1):
			if res in cat_res:
				# Exclude the catalytic residues from the designable list
				continue
			a_cords = res_ca_cords(pose, res)
			for i in pep_coords:
				# If any res is within radius of any pep res, add to the list
				if point_dist(a_cords, i) <= radius:
					mutable_residues.append(res)
					break

	return pep_start, pep_len, mutable_residues 


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


def quick_thread(pdb, sequences, destination):
	pose = pose_from_pdb(pdb)

	out_names = []

	for seq in sequences:
		threaded_pose = Pose()
		threaded_pose.assign(pose)
		threaded_pose = thread_seq(threaded_pose, 197, 11, seq)

		# Outputting model
		out_name = join(destination, pdb.replace('.pdb', '_' + seq + '.pdb.gz'))
		threaded_pose.dump_pdb(out_name)
		out_names.append(out_name)

	return out_names


def apply_constraints(pose):
	""" Applies the constraints form the input CST file to a pose """
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)
	return pose


def fastrelax(pose, near_res, score_function):
	""" Runs the FastRelax protocol on a pose, using the default score """
	relax = FastRelax()
	relax.set_scorefxn(score_function)

	# Movemap
	mm = MoveMap()
	mm.set_bb_true_range(197,207)
	mm.set_chi_true_range(197,207)
	for i in near_res:
		mm.set_chi(i, True)
	relax.set_movemap(mm)

	relax.apply(pose)
	return pose


def set_relax(pdb, num_decoys, des_res, score_function, dest):
	"""
	Outputs a set number of relaxed decoys, and returns the single lowest 
	scoring model
	"""
	pose = apply_constraints(pose_from_pdb(pdb))
	sf = score_function
	dec_name = pdb.replace('.pdb.gz', '_relaxed')
	jd = PyJobDistributor(dec_name, num_decoys, sf)

	best_structure = Pose()
	best_structure.assign(pose)
	best_score = sf(pose)

	while not jd.job_complete:
		pp = Pose()
		pp.assign(pose)
		pp = fastrelax(pp, des_res, sf)
		pose_score = sf(pp)

		if pose_score < best_score:
			best_score = pose_score
			best_structure.assign(pp)

		jd.output_decoy(pp)

	best_structure.dump_pdb(pdb.replace('.pdb', '_best.pdb'))
	return best_structure


def thread_and_relax(pdb, pep_sequences, dir_nam, des_res, sf):
	"""
	Threads given peptide sequences onto the starting structure. Then uses 
	the set_relax function to generate a set of relaxed structures for each,
	identifying the one with the lowest score as the best. Decoys are output
	into a specified subdirectory.
	"""
	# Threading the peptide sequences
	out_files = quick_thread(pdb, pep_sequences, dir_nam)

	# For each sequence, creating relax decoys and returning best one
	for i in out_files:
		# Producing ten relaxed decoys, storing the best one
		top_decoy = set_relax(i, 11, des_res, sf, dir_nam)


def get_best_structs(source_dir):
	"""
	Returns a list of best structures in an input directory (_best.pdb).
	PDB files are loaded separately, rather than being passed from a previous
	setp, because the constraints should not be applied for the ddG scoring.
	"""
	file_loc = join(basename(source_dir),"*best.pdb.gz")
	bests = glob(file_loc)
	best_structures = [pose_from_pdb(i) for i in bests]

	return best_structures


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

	# Packer for protease
	pro_task = standard_packer_task(pro_pose)
	pro_task.restrict_to_repacking()
	pro_task.temporarily_fix_everything()
	for i in near_res:
		pro_task.temporarily_set_pack_residue(i, True)

	# Packer for peptide
	pep_task = standard_packer_task(pep_pose)
	pep_task.restrict_to_repacking()

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


def preappend_lines(f_name, lines):
	""" Adds text at the beginning of a file """
	with open(f_name, 'r+') as f:
		orig = f.read()
		f.seek(0,0)
		f.write(lines + '\n' + orig)


def relaxed_score(pose, near_res, score_function):
	""" 
	Performs two fast relaxes, one on the protease, and the other on the 
	peptide. This allows for the determination og the ddG. 
	Movemap matches those used in the set relax for consistency.
	"""
	# Splitting pose
	chains = pose.split_by_chain()
	prot_pose = Pose()
	prot_pose.assign(chains[1])
	pep_pose = Pose()
	pep_pose.assign(chains[2])

	sf = score_function
	relax = FastRelax()
	relax.set_scorefxn(sf)

	# Movemap, protease
	mm_prot = MoveMap()
	mm_prot.set_bb(False)
	mm_prot.set_chi(False)
	for i in near_res:
		mm_prot.set_chi(i, True)

	# Movemap, peptide
	mm_pep = MoveMap()
	mm_pep.set_bb(True)
	mm_pep.set_chi(True)

	# Relaxing and scoring
	relax.set_movemap(mm_prot)
	relax.apply(prot_pose)
	prot_score = sf(prot_pose)

	relax.set_movemap(mm_pep)
	relax.apply(pep_pose)
	pep_score = sf(pep_pose)

	# Getting scores
	dock_score = sf(pose)
	split_score = prot_score + pep_score 
	ddg = dock_score - split_score

	return dock_score, split_score, ddg 


def add_const_remarks(r1_chain, r1_no, r1_name, r2_chain, r2_no, r2_name, cst_num):
	""" 
	Formats a remark to apply a constraint in a PDB file, given two residue
	chain numbers, residue unobers, and residue 3-letter names. 
	"""
	template = "REMARK 666 MATCH TEMPLATE {} {} {} \
				MATCH MOTIF {} {} {}  1  {}\n"
	remark_line = template.format(  r1_chain, str(r1_no), r1_name, 
									r2_chain, str(r2_no), r2_name,
									str(cst_num))
	return remark_line


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
	dir_nam = pdb.replace('.pdb', '_decoys')
	if not isdir(dir_nam):
		makedirs(dir_nam)

	# Determining peptide part of PDB file, residues near peptide
	pep_start, pep_length, des_res = \
		res_to_design(pdb, args.radius_design, args.cat_res)

	# Reading inputs for peptide sequences
	pep_sequences = get_seq_list(args.peptide_sequence)

	# Cresting structures, if this is a first use of the script
	if not args.rerun:
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


if __name__ == '__main__':
	main()