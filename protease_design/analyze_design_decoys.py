"""
Write a good docstring
"""

import argparse
import design_protease as dp
from glob import glob
import gzip
#from numpy import mean
#from operator import itemgetter as iget
from os.path import basename, join, isfile
import pickle as pickle
from pyrosetta import *
from pyrosetta.rosetta.core.simple_metrics.metrics import \
	InteractionEnergyMetric, RMSDMetric, SasaMetric, \
	SequenceSimilarityMetric, TotalEnergyMetric
from pyrosetta.rosetta.protocols.simple_ddg import ddG
#from pyrosetta.teaching import EMapVector
import subprocess

def parse_args():
	info = "Analyze design decoys from design_protease.py"
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--folder", type=str, required=True, 
		help="Pick folder to check")
	parser.add_argument("-des", "--design_type", type=str, required=True, 
		default='e', choices=['e', 's', 'b'], help="What kind of design was \
		used on this set? Options: 'e' for enzyme-only, 's' for \
		substrate-only, 'b' for both. Default is enzyme-only.")
	parser.add_argument("-cons", "--constraints", type=str, required=True, 
		default='ly104.cst', help="Pick constraints file")
	parser.add_argument("-hbn", "--use_hb_net", action="store_true", 
		help="Option to include HBnet score term in design.")
	parser.add_argument("-eval", "--force_evaluate", action = "store_true", 
		help="Re-evaluate decoys")
	parser.add_argument("-q", "--quiet", action="store_true", 
		help="Do you want to mute Rosetta output")
	args = parser.parse_args()
	return args


def unzip_file(fi):
	""" 
	Read in a gzipped file and save an unzipped version, stripping off .gz 
	from the file name ending, and returning that rename.
	"""
	# Unzip the file
	subprocess.call(['gunzip', fi])

	# Remove .gz from filename
	assert fi[-3:] == '.gz'
	out_file = fi[:-3]

	return out_file


def zip_file(fi):
	""" 
	Read in an unzipped file and save compress it, adding .gz to the file name 
	ending, and returning that rename.
	"""
	# Unzip the file
	subprocess.call(['gzip', fi])

	# Add .gz to filename
	assert fi[-3:] != '.gz'
	out_file = fi + '.gz'

	return out_file


def fix_file(pdb):
	"""
	Cleans up a PDB file, stripping out LINK lines and changing enzdes headers 
	so that constraints won't crash. Overwrites the original with the cleaned 
	version.
	"""
	with open(pdb, 'r') as r:
		lines = r.readlines()

	# Collect lists of LINK lines, enzdes header lines, and range of atom lines
	link_lines = []
	enzdes_lines = []
	first_atom = 0
	last_atom = 0
	for n, line in enumerate(lines):
		# Link lines
		if 'LINK' in line:
			link_lines.append(n)

		# Enzdes lined
		if 'REMARK 666' in line:
			enzdes_lines.append(n)

		# Atom lines range
		if 'ATOM' in line and first_atom == 0:
			first_atom = n
		if 'ATOM' in line and n > last_atom:
			last_atom = n

	# Remove LINK lines, which interfere with putting in constraints
	for l2r in link_lines:
		lines[l2r] = '\n'

	# Fixing enzdes comments block so comments match (mutated) PDB sequence
	for eline in enzdes_lines:
		line_text = lines[eline]
		# Splitting into columns. Columns 4-6 refer to the first residue, and 
		# columns 9-11 refer to the second residue
		e_columns = line_text.split()

		# Checking whether enzdes constraint headers match sequence in PDB
		for e_chain, e_aa, e_res in [e_columns[4:7], e_columns[9:12]]:
			orig_text = ' '.join([e_chain, e_aa, e_res])

			# Looping through all ATOM lines in the PDB until finding the 
			# right residue, then checking agreement and stopping loop
			for atom_line in lines[first_atom: last_atom]:
				# Skip any line that isn't an atom line, such as a TER
				if ('ATOM' not in atom_line):
					continue

				# Splitting into columns. Column 4 is the chain, column 5 is 
				# the residue number, column 3 is the aa name
				a_columns = atom_line.split()
				a_aa, a_chain, a_res = a_columns[3:6]

				# Skip until finding the right residue in the right chain
				if (a_chain != e_chain) or (a_res != e_res):
					continue

				# If the enzdes header and PDB agree, stop there
				if e_aa == a_aa:
					break

				# If the enzdes header and PDB disagree, correct the header
				else:
					correct_text = ' '.join([a_chain, a_aa, a_res])
					line_text = line_text.replace(orig_text, correct_text)
					break

		# Setting header line with corrected text
		lines[eline] = line_text
		
	# Saving corrected PDB
	with open(pdb,'w') as w:
		w.writelines(lines)


class single_design_analysis():
	"""
	Data storage object for a design decoy. Opens a design decoy and it's 
	relax (pre-design) counterpart and compares the two. Stores a list of 
	mutable residues and whether they did mutate, and what to. Tracks changes
	to the peptide sequence. Collects scoring information on the designed 
	decoy as well, including total energy, constraint energy, and shell energy, 
	and calculating a discriminator value.
	"""
	def __init__(self, relax_pdb, design_pdb, design_type):
		# PDB file names
		self.relax_pdb = relax_pdb
		self.design_pdb = design_pdb

		# Residue identifications (peptide, protease)
		self.peptide_residues = []
		self.original_peptide_sequence = ''
		self.final_peptide_sequence = ''

		self.mutable_residues_pose = []
		self.mutable_residues_pdb = []
		self.mutated_residues = {}

		# Energies
		self.total_energy = 0
		self.protease_total_energy = 0
		self.peptide_total_energy = 0
		self.interaction_energy = 0
		self.ddg_energy = 0
		self.hbnet_energy = 0
		self.constraint_energy = 0
		self.sasa_burial = 0

		# Energy changes (designed - relaxed)
		self.total_energy_change = 0
		self.protease_total_energy_change = 0
		self.peptide_total_energy_change = 0
		self.interaction_energy_change = 0
		self.ddg_energy_change = 0
		self.hbnet_energy_change = 0
		self.constraint_energy_change = 0
		self.sasa_burial_change = 0
		self.peptide_sequence_similarity = 0
		self.protease_sequence_similarity = 0

		# Run comparison analysis
		self.compare_design(design_type)

	def make_selectors(self, design_type):
		""" 
		Makes residue selectors based on design_protease.py, determining which
		residues were designable similar to how the design protocol did.

		Catalytic residues are not excluded from the designable set, but 
		will not have changed. Full substrate peptide will be treated as 
		peptide, so any subset will be ignored.
		"""
		# Set appropriate design switches for selection
		if design_type in ['e', 'b']:
			des_protease = True
		else:
			des_protease = False

		if design_type in ['s', 'b']:
			des_peptide = True
		else:
			des_peptide = False

		# Create selectors
		selectors = dp.select_residues(None, None, 
			design_protease=des_protease, design_peptide=des_peptide)

		return selectors

	def calc_energies(self, designed_pose, relaxed_pose, selectors):
		"""
		Calculate energies for the designed pose, and changes in those energies
		between the designed pose and the relaxed pose. 
		"""
		# Total energies
		sf_default = dp.get_score_function(constraints=False)
		tem = TotalEnergyMetric()
		tem.set_scorefunction(sf_default)
		self.total_energy = tem.calculate(designed_pose)
		self.total_energy_change = \
			self.total_energy - tem.calculate(relaxed_pose)

		# Protease energies
		tem.set_residue_selector(selectors['protease'])
		self.protease_total_energy = tem.calculate(designed_pose)
		self.protease_total_energy_change = \
			self.protease_total_energy - tem.calculate(relaxed_pose)

		# Peptide energies
		tem.set_residue_selector(selectors['peptide'])
		self.peptide_total_energy = tem.calculate(designed_pose)
		self.peptide_total_energy_change = \
			self.peptide_total_energy - tem.calculate(relaxed_pose)

		# Interaction energies
		iem = InteractionEnergyMetric(selectors['protease'], 
										selectors['peptide'])
		self.interaction_energy = iem.calculate(designed_pose)
		self.interaction_energy_change = \
			self.interaction_energy - iem.calculate(relaxed_pose)

		# ddG
		ddg = ddG(sf_default, 1) # 1 specifies the jump number
		ddg.calculate(designed_pose)
		self.ddg_energy = ddg.sum_ddG()
		ddg.calculate(relaxed_pose)
		self.ddg_energy_change = self.ddg_energy - ddg.sum_ddG()

		# HBnet energies
		sf_hb = dp.get_score_function(
			ref15=False, constraints=False, hbnet=True)
		#tem_hb = TotalEnergyMetric()
		#tem_hb.set_scorefunction(sf_hb)
		#self.hbnet_energy = tem_hb.calculate(designed_pose)
		#self.hbnet_energy_change = \
		#	self.hbnet_energy - tem_hb.calculate(relaxed_pose)
		self.hbnet_energy = sf_hb(designed_pose)
		self.hbnet_energy_change = sf_hb(relaxed_pose)

		# Constraint energies
		sf_cst = dp.get_score_function(ref15=False)
		tem_cst = TotalEnergyMetric()
		tem_cst.set_scorefunction(sf_cst)
		des_pose_cst = Pose(designed_pose)
		dp.apply_constraints(des_pose_cst)
		rel_pose_cst = Pose(relaxed_pose)
		dp.apply_constraints(rel_pose_cst)
		self.constraint_energy = tem_cst.calculate(des_pose_cst)
		self.constraint_energy_change = \
			self.constraint_energy - tem_cst.calculate(rel_pose_cst)

		# SASA burial
		sasa = SasaMetric()
		sasa.set_residue_selector(selectors['peptide'])
		d_docked_sasa = sasa.calculate(designed_pose)
		d_undocked_sasa = sasa.calculate(designed_pose.split_by_chain()[2])
		self.sasa_burial = d_undocked_sasa - d_docked_sasa
		r_docked_sasa = sasa.calculate(relaxed_pose)
		r_undocked_sasa = sasa.calculate(relaxed_pose.split_by_chain()[2])
		r_sasa_burial = r_undocked_sasa - r_docked_sasa
		self.sasa_burial_change = self.sasa_burial - r_sasa_burial

		# Sequence similarities
		ssm = SequenceSimilarityMetric()
		ssm.set_native_pose(relaxed_pose)
		ssm.set_residue_selector(selectors['peptide'])
		self.peptide_sequence_similarity = ssm.calculate(designed_pose)
		ssm.set_residue_selector(selectors['mutable'])		
		self.protease_sequence_similarity = ssm.calculate(designed_pose)


	def compare_design(self, design_type):
		"""
		"""
		# Make selectors
		selectors = self.make_selectors(design_type)

		# Open PDBs
		designed_pose = pose_from_pdb(self.design_pdb)
		relaxed_pose = pose_from_pdb(self.relax_pdb)

		# Get peptide residues list and peptide sequences
		self.peptide_residues = dp.selector_to_list(relaxed_pose, 
			selectors['peptide'])
		for sr in self.peptide_residues:
			self.original_peptide_sequence += relaxed_pose.residue(sr).name1()
			self.final_peptide_sequence += designed_pose.residue(sr).name1()

		# Get designable residues list (both pose and PDB)
		self.mutable_residues_pose = dp.selector_to_list(relaxed_pose, 
			selectors['mutable'])
		rlx_inf = relaxed_pose.pdb_info()
		mrp = [rlx_inf.pose2pdb(er) for er in self.mutable_residues_pose]
		# pose2pdb outputs a string with the res number then the chain letter
		self.mutable_residues_pdb = [int(i.split()[0]) for i in mrp]

		# Identify changes in protease sequence
		for n, i in enumerate(self.mutable_residues_pose):
			orig_aa = relaxed_pose.residue(i).name1()
			final_aa = designed_pose.residue(i).name1()

			if final_aa == orig_aa:
				self.mutated_residues[self.mutable_residues_pdb[n]] = \
					[orig_aa, 'unchanged']

			else:
				self.mutated_residues[self.mutable_residues_pdb[n]] = \
					[orig_aa, final_aa]

		# Get scoring information
		self.calc_energies(designed_pose, relaxed_pose, selectors)


def analyze_design_decoys(args):
	"""
	Reading in all the relaxed/designed PDBs in a folder, this function 
	collects various pieces of information from each, storing that information 
	in a list of class objects.
	"""
	# Get list of designed pdbs (even if gzipped)
	designed_pdbs = glob(join(args.folder, '*designed*.pdb*'))
	designed_pdbs.sort()

	# Loop through all PDBs, unzipping if necessary, skipping incomplete files
	decoys_data_extract = []
	siz = len(designed_pdbs)
	skip_count = 0
	for n, d in enumerate(designed_pdbs):
		print('\n', "{} / {}\t{}".format(n, siz, basename(d)))
		# Skip incomplete 'in_progress' files
		if 'in_progress' in d:
			print("\nWarning: skipped {}\n".format(d))
			skip_count += 1
			continue

		# Identifying relaxed version
		dpdb = d
		rpdb = dpdb.replace('_designed_', '_relaxed_')

		# Unzip if gzipped
		if '.gz' in dpdb:
			dpdb = unzip_file(dpdb)
			rpdb = unzip_file(rpdb)

		# Cleaning up PDB files
		fix_file(dpdb)
		fix_file(rpdb)

		# Create data object for this decoy, add it to the list
		ddata = single_design_analysis(rpdb, dpdb, args.design_type)
		decoys_data_extract.append(ddata)

		# Rezip files
		dpdb = zip_file(dpdb)
		rpdb = zip_file(rpdb)

	print("Skipped {} in_progress files".format(skip_count))
	return decoys_data_extract


def write_summary(filename, design_analysis_list, peptide=False):
	"""
	Write summary csv file
	"""
	header = ['name', 'substrate_sequence', 'total_energy', 'protease_energy',
		'peptide_energy', 'constraints', 'hbnet', 'ddg', 'interaction_energy', 
		'sasa_burial', 'peptide_similarity', 'protease_similarity',
		'total_energy_change', 'protease_change', 'peptide_change', 
		'constraints_change', 'hbnet_change', 'ddg_change', 
		'interaction_change', 'sasa_burial_change']

	# Make mutations list header
	mutated_sites = []
	site_origs = []
	for p in design_analysis_list:
		for k, v in p.mutated_residues.items():
			if k not in mutated_sites:
				mutated_sites.append(k)
				site_origs.append(v[0])

	site_origs = [x for _,x in sorted(zip(mutated_sites,site_origs))]
	mutated_sites.sort()
	for n, i in enumerate(mutated_sites):
		header.append(site_origs[n] + str(i))


	for p in design_analysis_list:
		if p == design_analysis_list[0]: # First decoy
			with open(filename, 'w') as o:
				str_head = ','.join(header) + '\n'
				o.write(str_head)

		line = []
		line.append(basename(p.design_pdb))
		line.append(p.final_peptide_sequence)
		line.append(p.total_energy)
		line.append(p.protease_total_energy)
		line.append(p.peptide_total_energy)
		line.append(p.constraint_energy)
		line.append(p.hbnet_energy)
		line.append(p.ddg_energy)
		line.append(p.interaction_energy)
		line.append(p.sasa_burial)
		line.append(p.peptide_sequence_similarity)
		line.append(p.protease_sequence_similarity)
		line.append(p.total_energy_change)
		line.append(p.protease_total_energy_change)
		line.append(p.peptide_total_energy_change)
		line.append(p.constraint_energy_change)
		line.append(p.hbnet_energy_change)
		line.append(p.ddg_energy_change)
		line.append(p.interaction_energy_change)
		line.append(p.sasa_burial_change)

		for ms in mutated_sites:
			if ms in p.mutated_residues:
				line.append(p.mutated_residues[ms][1][0])
			else:
				line.append('')	

		with open(filename, 'a') as o:
			str_line = ','.join([str(i) for i in line]) + '\n'
			o.write(str_line)

	print('Saved as {}'.format(filename))

def main(args):
	# Determining target folder with decoys to analyze
	folder = args.folder
	folder_base = basename(folder.rstrip('/'))

	# Naming pickled analysis results
	pickle_name = join(folder, folder_base + '_design_analysis.pkl')

	# Reading in previous analysis if it has been done, unless user said redo
	if isfile(pickle_name) and not args.force_evaluate:
		# Retrieving pickled results
		with open(pickle_name, 'rb') as i:
			collected_decoy_info = pickle.load(i)

	# If analysis needs to be performed, perform it and pickle the results
	else:
		# Performing analysis
		print("\n\nAnalyzing\n\n")
		collected_decoy_info = analyze_design_decoys(args)

		# Pickling results
		with open(pickle_name, 'wb') as o:
			pickle.dump(collected_decoy_info, o, protocol=-1)
		print('\nRaw results saved as {}'.format(pickle_name))

	# Write report
	report_name = join(folder, folder_base + '_summary.csv')
	write_summary(report_name, collected_decoy_info)


if __name__ == '__main__':
	# Parse user inputs
	args = parse_args()

	# Initialize PyRosetta
	ros_opts = '-enzdes::cstfile {} '.format(args.constraints)
	ros_opts += '-run:preserve_header '
	if args.quiet:
		ros_opts += '-mute all '
	init(ros_opts)

	# Execute script
	main(args)


'''
Bugging on selectors--non-uniform
Pick unique decoys, compare decoys with BLOSSUM
#Fix hbnet scoring
#SASA of pocket residues
Interface H-bonds count
#Seq difference between design and original
'''