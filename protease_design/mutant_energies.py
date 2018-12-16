#!/usr/bin/python
""" 
PyRosetta4, Python 2.7
Joseph Lubin, 2017
"""
import argparse
from glob import glob
from numpy import mean
from operator import itemgetter as iget
from os.path import basename, join, isfile
from pyrosetta import *
from pyrosetta.teaching import EMapVector
from pyrosetta.rosetta import *
from design_protease import *
import cPickle as pickle

# Residue type positive (+), negative (-), neutral polar (N), hydrophobic (O)
res_types = {'A':'O', 'C':'N', 'D':'-', 'E':'-', 'F':'O', 'G':'O', 
			'H':'+', 'I':'O', 'K':'+', 'L':'O', 'M':'O', 'N':'N', 'P':'O', 
			'Q':'N', 'R':'+', 'S':'N', 'T':'N', 'V':'O', 'W':'O', 'Y':'O'}

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("folder", help="Pick folder to check")
	parser.add_argument("-des_pep", "--design_peptide", action = "store_true", 
		help="Expand residue selection to peptide, not just protease")
	parser.add_argument("-eval", "--force_evaluate", action = "store_true", 
		help="Re-evaluate decoys")
	args = parser.parse_args()
	return args


class mutation_collection:
	"""
	This class object will include information on the mutations in decoys 
	produced by design_protease.py
	"""
	def __init__(self, folder, sequence, design_peptide=False):
		# Identifying decoy location
		self.folder = folder

		# Identifying peptide sequence, and whether it is cleaved
		self.sequence = sequence
		self.seq_and_cleavage()

		# Getting lists of designed and relaxed decoys
		self.collect_models()

		# Defining relevant residues
		self.design_peptide = design_peptide
		self.select_residues()

		# Setting score function
		self.sf = get_fa_scorefxn()

		# Identifying mutations and their frequency and energies
		other_class_values = ['mutations', 'mut_freq', 'mutation_pdb', 
			'mut_location', 'mut_res_es', 'mut_int_es', 
			'interaction_residues']
		for i in other_class_values:
			setattr(self, i, [])

		# Iterating through all decoys extracting relevant info
		for i in range(self.decoy_count):
			print '\tanalyzing decoy', i + 1
			self.ident_mutations(self.relaxed_pdbs[i], self.designed_pdbs[i])

		# Getting mutation frequencies, energy averages
		self.mut_rate = \
			[round(float(i)/self.decoy_count,3) for i in self.mut_freq]
		self.mut_res_min_e = [min(i) for i in self.mut_res_es]
		self.mut_res_average_e = [mean(i) for i in self.mut_res_es]
		self.mut_int_min_e = [min(i) for i in self.mut_int_es]
		self.mut_int_average_e = [mean(i) for i in self.mut_int_es]

	##########################################################################
	def classify_residues(self):
		""" 
		Looks at the residues in a sequence and determines whether they are 
		positively charged (+), negatively charged (-), neutral-polar (N), or
		hydrophobic (O).
		"""
		self.peptide_res_types = ''
		self.peptide_charge = 0
		self.peptide_polarity = 0

		for res in self.short_sequence:
			res_type = res_types[res]
			self.peptide_res_types += res_type

			if res_type == '+':
				self.peptide_charge += 1
			if res_type == '-':
				self.peptide_charge -= 1

			if res_type in ['+', '-', 'N']:
				self.peptide_polarity += 1

	def seq_and_cleavage(self):
		""" 
		Determines the input peptide sequence and whether it was from a 
		cleaved or uncleaved substrate. The sequences fed into this object
		by isolate_sequence_from_fasc will have either SEQUENCE or unSEQUENCE,
		depending on whether the sequence is cleaved or not.
		"""
		if 'un' in self.sequence:
			self.sequence = self.sequence.replace('un','')
			self.cleaved = 'uncleaved'
		else:
			self.cleaved = 'cleaved'
		print '\n', self.sequence, self.cleaved

		self.short_sequence = self.sequence[1:7]
		self.classify_residues()
		print self.short_sequence

	def collect_models(self):
		""" Collects threaded, designed, and relaxed decoys """
		# Finding the threaded pre-relax model
		thread_search = join(self.folder, '*' + self.sequence + '.pdb.gz')
		threaded_structures = glob(thread_search)
		assert len(threaded_structures) == 1 
			# There should only be one threaded model for each sequence
		self.threaded_pdb = threaded_structures[0]

		# Finding relaxed decoys
		relaxed_search = join(self.folder, '*' + self.sequence + '_relaxed*.pdb*')
		self.relaxed_pdbs = sorted(glob(relaxed_search))

		# Finding designed decoys
		designed_search = join(self.folder, '*' + self.sequence + '_designed*.pdb*')
		self.designed_pdbs = sorted(glob(designed_search))

		# Counting decoys
		assert len(self.designed_pdbs) == len(self.relaxed_pdbs)
		self.decoy_count = len(self.designed_pdbs)

	def select_residues(self):
		""" 
		Uses residue selectors to determine the set of designable residues.
		mutable_residues_selector and selector_to_list are from 
		design_protease. The option to select peptide residues is passed as 
		true so that any mutations will be collected, regardless of whether 
		the peptide was designable.
		"""
		peptide_selector = ResidueIndexSelector('198-203')
		designable_selector = mutable_residues_selector(self.design_peptide)

		threaded_pose = pose_from_pdb(self.threaded_pdb)
		self.pep_res = selector_to_list(threaded_pose, peptide_selector)
		self.des_res = selector_to_list(threaded_pose, designable_selector)

	##########################################################################
	def residue_before_and_after(self, res, before, after):
		""" Determines the AA in the same position in two poses """
		start_res = before.residue(res)
		start_res_name = start_res.name1()

		end_res = after.residue(res)
		end_res_name = end_res.name1()

		mut_string = '_'.join([str(res), start_res_name, end_res_name])

		return start_res_name, end_res_name, mut_string

	def res_interact_energy(self, pose, res1, res2):
		""" 
		Takes a given pose and residue pair, and sums the weighted residue-pair 
		interaction energies between the two.
		ref_wts from design_protease.
		"""
		# Collecting interaction energies
		sum_int_energy = 0
		r1 = pose.residue(res1)
		r2 = pose.residue(res2)

		# Getting emap
		emap = EMapVector()
		self.sf.eval_ci_2b(r1, r2, pose, emap) # Includes many zero terms

		# Cleaning up emap
		org_emp = {}
		emap_nz = str(emap.show_nonzero()).split()
		for i in range(len(emap_nz)/2): # Converting to a dict
			org_emp[emap_nz[2 * i].rstrip(':')] = float(emap_nz[2 * i + 1])

		# Adding weighted scores
		for sc_term in org_emp:
			if sc_term in ref_wts: # Emap contains terms with 0 weight in REF
				sum_int_energy += org_emp[sc_term] * ref_wts[sc_term]

		return sum_int_energy

	def residue_pairwise_energy(self, res, before, after):
		""" 
		Theis function determines whether a mutation is on the protease or  
		peptide, and conversely the relevant set of residues contributing to 
		the mutated residue's score. Then it runs interact_energy accordingly
		on the residue in the relaxed and designed poses, and calculates the
		difference. Returns the location of the mutation, the residue 
		interaction energy change between the relaxed and designed models, 
		and a list of the residues for which the interaction energy was 
		nonzero.
		"""
		if res in self.pep_res:
			mut_loc = 'Peptide'
			partner = self.des_res
		else:
			mut_loc = 'Protease'
			partner = self.pep_res

		# Getting sum residue interaction energy with contact residues
		# and residue contact list
		int_e_sum = 0
		contacts = {}
		for n, p in enumerate(partner):
			res_start_int_e = self.res_interact_energy(before, res, p)
			res_end_int_e = self.res_interact_energy(after, res, p)
			res_delta_int_e = res_end_int_e - res_start_int_e
			int_e_sum += res_delta_int_e
			contacts[p] = [res_start_int_e, res_end_int_e]

		return mut_loc, int_e_sum, contacts

	def ident_mutations(self, rel_pdb, des_pdb):
		"""
		Compares the sequences of a decoy before and after design at specified 
		residues, and identifies the differences. Mutations are identified in 
		the format of N_A_B, where N is the residue number, A is the starting 
		residue, and B is the ending residue. The residue energy difference for 
		the mutation is also collected, as is the sum of two-residue 
		interactions between the mutated residue and the peptide residues (if 
		the mutated residue is on the protease) or between the mutated 
		residue and the mutable protease residues (if the residue is on the 
		peptide). res_scores ia a function from design_protease.
		"""
		# Getting starting and end poses, excerpts from residue energy tables
		relax_pose = pose_from_pdb(rel_pdb)
		design_pose = pose_from_pdb(des_pdb)
		self.sf(relax_pose)
		self.sf(design_pose)

		# Evaluating protease mutations
		for i in self.des_res:
			# Selecting residues to compare
			start_res_name, end_res_name, mut_string = \
				self.residue_before_and_after(i, relax_pose, design_pose)
			
			# If they are different, i.e. there was a mutation
			if start_res_name != end_res_name:
				# Getting overall residue energy change
				relax_res_energy = res_scores(relax_pose, [i], self.sf)[0]
				design_res_energy = res_scores(design_pose, [i], self.sf)[0]
				res_delta_e = design_res_energy - relax_res_energy

				# Getting pairwise residue energy change
				mut_loc, res_delta_int_e, interact_res = \
					self.residue_pairwise_energy(i, relax_pose, design_pose)

				# Aggregating results
				if mut_string in self.mutations:
					# Mutation has been identified on previous decoy
					mut_ind = self.mutations.index(mut_string)
					assert self.mut_location[mut_ind] == mut_loc
					self.mut_freq[mut_ind] += 1
					self.mutation_pdb[mut_ind].append(des_pdb)
					self.mut_res_es[mut_ind].append(res_delta_e)
					self.mut_int_es[mut_ind].append(res_delta_int_e)
					self.interaction_residues[mut_ind].append(interact_res)
				else:
					# New mutation or first decoy
					self.mutations.append(mut_string)
					self.mut_location.append(mut_loc)
					self.mut_freq.append(1)
					self.mutation_pdb.append([des_pdb])
					self.mut_res_es.append([res_delta_e])
					self.mut_int_es.append([res_delta_int_e])
					self.interaction_residues.append([interact_res])


def isolate_sequence_from_fasc(fasc_name):
	""" 
	Reads the name of a FASC file output by design_protease.py, isolates and 
	returns the peptide sequence 
	"""
	short_name = basename(fasc_name)
	minus_tail = short_name.replace('_designed.fasc', '')
	sequence = minus_tail.replace('cleaved_ly104_wt_','')
	return sequence


def mutations_by_seq_report(name, mc, head=False):
	"""
	Writes a mutation_summary report, where the mutations for each decoy set
	are summarized.
	"""
	# Formatting report content
	template = '{:10s}' + '{:15s}' + '{:10s}' * 3 + '{:18s}' * 2 + '\n'
	
	# Writing header
	if head:
		header = ['Cleaved', 'Sequence', 'Location', 'Mutation', 'Frequency', 
				'delta_E_residue', 'delta_E_interaction']
		with open(name, 'w') as r:
			r.write(template.format(*header))

	# Writing body
	report_lines = []
	for i in range(len(mc.mutations)):
		line = [str(j) for j in [mc.cleaved, mc.short_sequence, 
				mc.mut_location[i], mc.mutations[i], mc.mut_rate[i], 
				round(mc.mut_res_min_e[i], 3), 
				round(mc.mut_int_min_e[i], 3)]]
		report_lines.append(line)
	report_lines.sort()

	with open(name, 'a') as r:
		for l in report_lines:
			r.write(template.format(*l))


class mutations_aggregate:
	"""
	This class object collects mutations information across a set of 
	mutation_collection objects and performs some analysis.
	"""
	def __init__(self, mc_set):
		self.mc_set = mc_set
		self.design_peptide = self.mc_set[0].design_peptide

		# Get full list of mutable residues (differs within the set), pep res
		self.residues_list = self.get_full_residue_list()
		self.mutable_count = len(self.residues_list)
		self.pep_res = self.get_pep_res()

		# Initializing data bins
		self.total_potential_models = 0
		self.actual_potential_models = [0] * self.mutable_count
		self.mutations_counts = [0] * self.mutable_count

		self.mutation_occurrences = [{} for i in range(self.mutable_count)] 
		self.mutation_res_total_energies_raw = \
			[{} for i in range(self.mutable_count)]
		self.mutation_res_interact_energies_raw = \
			[{} for i in range(self.mutable_count)]

		# Populating bins
		self.starting_residues = self.get_initial_residues()
		self.collect_mutations_and_frequencies()
		self.count_mutations()
		self.condense_energies()

		# Identifying residue interactions
		self.residue_interactions = self.get_residue_interactions()

		# Creating report file
		self.output_length = 0
		self.pre_report = self.mc_tabulation()
		self.report = self.cleanup_report_table()

		# Creating residue pairing representative decoys dict list
		self.representative_deocys = self.pick_representative_deocys()

	def get_full_residue_list(self):
		""" 
		Reads through a set of mutation_collection objects and collects a full
		list of all residues that were mutable across the simulations.
		""" 
		complete_residues_list = []
		for mc in self.mc_set:
			for i in mc.des_res:
				if i not in complete_residues_list:
					complete_residues_list.append(i)

		complete_residues_list.sort()
		return complete_residues_list

	def get_pep_res(self):
		""" 
		Checks through the set of mutation_collection objects, verifies that
		all have the same peptide residues identified, returns the set if so.
		"""
		pep_res = None

		for mc in self.mc_set:
			pep_list = mc.pep_res
			if pep_res == None: # Initial case
				pep_res = pep_list
			else:
				assert pep_res == pep_list

		return pep_res

	def get_initial_residues(self):
		""" Collect starting residues, before mutation. Sloppy hard-coded. """
		pose_seq = 'GSVVIVGRIILSGRGGPITAYAQQTRGLLGCIITSLTGRDKNQVEGEVQIVSTAAQT'
		pose_seq += 'FLATCINGVCWTVYHGAGTRTIASPKGPVIQMYTNVDQDLVGWPASQGTRSLTPCT'
		pose_seq += 'CGSSDLYLVTRHADVIPVRRRGDSRGSLLSPRPISYLKGSSGGPLLCPAGHAVGIF'
		pose_seq += 'RAAVCTRGVAKAVDFIPVENLETTMRS'

		starting_residues = []
		for i in self.residues_list:
			if i < len(pose_seq):
				starting_residues.append(pose_seq[i-1])

		# Matching length
		while len(starting_residues) < len(self.residues_list):
			starting_residues.append('var')

		return starting_residues

	def mc_name_to_table(self, mc):
		"""
		mutation_collection objects store mutations in the format of N_A_B, 
		where N is the residue number, A is the starting residue, and B is 
		the ending residue. For this application, a list of lists is more 
		useful. This function makes that conversion and appends the frequency 
		of mutation. Returns a list in form [N, B, F, E, I], where F is the 
		frequency with which the mutation was observed, E is the set of 
		total residue energy changes associated with the mutation, anf I is 
		the set of residue interaction energy changes associated with the 
		mutation.
		"""
		converted_mutations_list = [i.split('_') for i in mc.mutations]

		for n, i in enumerate(converted_mutations_list):
			 i[0] = int(i[0])
			 i.pop(1) # Remove original residue
			 i[1] = str(i[1])
			 i.append(mc.mut_freq[n])
			 i.append(mc.mut_res_es[n])
			 i.append(mc.mut_int_es[n])

		converted_mutations_list.sort()
		return converted_mutations_list

	def collect_mutations_and_frequencies(self):
		""" 
		Reads through the list of mutation_collection objects to extract a 
		list of mutations at each position, and also the total number of 
		decoys both within the set and which have each mutation.
		"""
		for mc in self.mc_set:
			# Getting sum count of all models in all sets
			self.total_potential_models += mc.decoy_count

			# Getting count of models with the potential for mutations, since
			# some model selectors allowed different residues
			for n, i in enumerate(self.residues_list):
				if i in mc.des_res:
					self.actual_potential_models[n] += mc.decoy_count

			# Collecting mutations
			mc_muts = self.mc_name_to_table(mc)
			for i in mc_muts:
				ind = self.residues_list.index(i[0])

				if i[1] in self.mutation_occurrences[ind]:
					self.mutation_occurrences[ind][i[1]] += i[2]
					self.mutation_res_total_energies_raw[ind][i[1]] += i[3]
					self.mutation_res_interact_energies_raw[ind][i[1]] += i[4]
				else:
					self.mutation_occurrences[ind][i[1]] = i[2]
					self.mutation_res_total_energies_raw[ind][i[1]] = i[3]
					self.mutation_res_interact_energies_raw[ind][i[1]] = i[4]

	def count_mutations(self):
		""" 
		Gets count of how many times a residue has mutated across all sets.
		Also makes list of rates of mutation occurrence 
		"""
		for n, i in enumerate(self.mutation_occurrences):
			for m in i:
				self.mutations_counts[n] += i[m]

		self.res_mutation_propensity = []
		for i in range(self.mutable_count):
			numerator = float(self.mutations_counts[i])
			denominator = float(self.actual_potential_models[i])
			propensity = round(numerator / denominator, 3)
			self.res_mutation_propensity.append(propensity)

	def condense_energies(self):
		""" 
		Lists of energies are associated with each mutation. This function
		gets the average and minimum of those sets.
		"""
		bin_set = ['mutation_res_tot_E_average', 'mutation_res_int_E_average',
					'mutation_res_tot_E_min', 'mutation_res_int_E_min']
		for i in bin_set:
			bin = [{} for x in range(self.mutable_count)]
			setattr(self, i, bin)

		for n, pos in enumerate(self.mutation_occurrences):
			for m in pos:
				self.mutation_res_tot_E_average[n][m] = \
					mean(self.mutation_res_total_energies_raw[n][m])
				self.mutation_res_tot_E_min[n][m] = \
					min(self.mutation_res_total_energies_raw[n][m])
				self.mutation_res_int_E_average[n][m] = \
					mean(self.mutation_res_interact_energies_raw[n][m])
				self.mutation_res_int_E_average[n][m] = \
					min(self.mutation_res_interact_energies_raw[n][m])

	def get_residue_interactions(self):
		""" 
		Looks at each mutation_collection and runs a residue 
		selector specific to each residue on the peptide, then converts 
		from a list of what mutable residues interact with each peptide 
		residue to a list of with which peptide residues each mutable residue 
		can interact. mutable_residues_selector and selector_to_list are from 
		design_protease.
		"""
		peptide_residues = range(198, 204)
		res_selectors = [mutable_residues_selector(single_pep_res=i) 
			for i in peptide_residues]
		pep_res_partners = [[] for i in peptide_residues]

		# Looking at each mutation_collection and running a residue selector
		for mc in self.mc_set:
			pose = pose_from_pdb(mc.threaded_pdb)
			for n, rs in enumerate(res_selectors):
				selection = selector_to_list(pose, rs)
				for i in selection:
					if i not in pep_res_partners[n]:
						pep_res_partners[n].append(i)

		# Switching from pep->mutable to mutable->pep
		int_size = self.mutable_count
		if self.design_peptide:
			int_size -= len(peptide_residues)
		residue_interactions = [[] for i in range(int_size)]

		for n, i in enumerate(pep_res_partners):
			i.sort()
			for j in i:
				ind = self.residues_list.index(j)
				residue_interactions[ind].append(peptide_residues[n])

		if self.design_peptide: 
			print residue_interactions + pep_res_partners
			return residue_interactions + pep_res_partners	
		else: 
			return residue_interactions	

	def string_interactions_together(self, position):
		""" Concatenates a list of interactable residues to a string """
		number_position = {'198':'p6', '199':'p5', '200':'p4',  
							'201':'p3', '202':'p2', '203':'p1'}

		partners = self.residue_interactions[position]
		partners_as_str = [str(i) for i in partners]

		for n, i in enumerate(partners_as_str):
			if i in number_position:
				partners_as_str[n] = number_position[i]

		partners_as_str.sort()

		return ','.join(partners_as_str)

	def back_convert_N_A_B(self, index, end_res):
		""" 
		Function to convert back to the mutation format in mutation_collection 
		"""
		locus = str(self.residues_list[index])
		initial_res = self.starting_residues[index]

		mut_string = '_'.join([locus, initial_res, end_res])
		return mut_string

	def extract_start_res_info(self, res, n_mut, n_tot):
		"""
		Pulls information about a residue/position prior to design
		"""
		mut_line_out = []
		mut_line_out.append(res)
		mut_line_out.append(res_types[res])
		occurrence = round(float(n_mut) / n_tot, 3)
		mut_line_out.append(str(occurrence))

		return mut_line_out

	def extract_end_mut_info(self, mutation_collection, mutation):
		""" 
		Pulls report information about a specific mutation from a 
		mutation_collection object.
		"""
		m_ind = mutation_collection.mutations.index(mutation)

		mc_extract = []
		mc_extract.append(mutation_collection.cleaved[0].upper())
		mc_extract.append(mutation_collection.short_sequence)
		t = '\'' + mutation_collection.peptide_res_types + '\''
		mc_extract.append(t)
		res_dE = round(mutation_collection.mut_res_min_e[m_ind], 3)
		mc_extract.append(str(res_dE))
		res_int_dE = round(mutation_collection.mut_int_min_e[m_ind], 3)
		mc_extract.append(str(res_int_dE))
		mc_extract.append(str(mutation_collection.mut_rate[m_ind]))

		return mc_extract

	def get_line_length(self, line):
		"""
		Determines the number of elements that will be written to a report 
		line. Each line is a list, though the last element may be a list of 
		lists. There are three layers. The output result is the number of 
		non-list elements across all three layers
		"""
		outer = sum([1 for i in line if type(i) is not list])
		middle = sum([1 for i in line[-1][0] if type(i) is not list])
		inner = len(line[-1][0][-1][0])
		return outer + middle + inner

	def mc_tabulation(self):
		""" 
		Loops through all identified mutations and for each, checks the list
		of mutation_collection objects for which ones yielded that mutation.
		For each collection that does include the mutation, this function 
		extracts some information about that collection's substrate sequence, 
		and about the statistics relating to the mutation within the 
		collection. Returns a list that gets processed for writing to a report 
		file.
		"""
		pre_output_list = []
		is_first = True
		for n, i in enumerate(self.mutation_occurrences):
			loc = self.residues_list[n]
			initial = self.starting_residues[n]
			propensity = self.res_mutation_propensity[n]
			potential = self.actual_potential_models[n]
			dr_for_all_models = (potential != self.total_potential_models)

			line_out = []
			line_out.append(str(loc))
			line_out.append(self.string_interactions_together(n))
			line_out.append(initial)

			if propensity == 0:
				# No mutations at this locus
				line_out.append('NO MUTATIONS'+ dr_for_all_models * '*')
				pre_output_list.append(line_out)

			else:
				if initial == 'var':
					line_out.append('var')
				else:
					line_out.append(res_types[initial])
				line_out.append(str(propensity) + dr_for_all_models * '*')
				line_out.append([])
				for j in i:
					mut_line_out = \
						self.extract_start_res_info(j, i[j], potential)
					mut_line_out.append([])	

					mc_mut = self.back_convert_N_A_B(n, j)
					for mc in self.mc_set:
						if mc_mut in mc.mutations:
							mc_extract = self.extract_end_mut_info(mc, mc_mut)
							mut_line_out[-1].append(mc_extract)

					line_out[-1].append(mut_line_out)

				pre_output_list.append(line_out)

				# Determining the length of the output string the first time through
				if is_first:
					self.output_length = self.get_line_length(line_out[:])
					is_first = False

		return pre_output_list

	def cleanup_report_table(self):
		""" 
		Converts table from mc_tabulation to a more readable form, intended
		to be opened in Excel, which leaves a cell visually empty for ="".
		"""
		cleaned_table = []
		for line in self.pre_report:
			# No mutations
			if any(['NO MUTATIONS' in i for i in line]):
				cleaned_table.append(line)

			# Mutations
			else:
				# Sorting
				line[-1] = sorted(line[-1], key=iget(1, 0), reverse=True)

				# Converting from list to lines
				for mutation in line[-1]:
					# Sorting
					mutation[-1] = sorted(mutation[-1], key=iget(0, 2))

					# Converting from list to lines
					for s in mutation[-1]:
						first_part = ['=""'] * len(line[:-1])
						second_part = ['=""'] * len(mutation[:-1])
						if s == mutation[-1][0]:
							second_part = mutation[:-1]
							if mutation == line[-1][0]:
								first_part = line[:-1]
						third_part = s

						cleaned_table.append(first_part + second_part + 
												third_part)

		return cleaned_table

	def favorable_interactions_by_pair_energy(self):
		"""
		Looks at the list of pairwise residue interactions across all 
		mutation_collection objects in the set, and collects a dict of decoys
		with protease mutations interacting with each peptide residue with an
		energy better than -0.5 KCal/mol. 
		"""
		interaction_sets = [{} for i in self.pep_res]

		for x, pep in enumerate(self.pep_res):
			for mc in self.mc_set:
				for y, mut in enumerate(mc.mutations):
					for z, en in enumerate(mc.interaction_residues[y]):
						e_start_round = round(en[pep][0],3)
						e_end_round = round(en[pep][1],3)
						e_dif = e_end_round - e_start_round
						#if e_end_round < -0.5:
						decoy_energy = [mc.short_sequence, e_end_round, 
										e_dif, mc.mutation_pdb[y][z]]
						if mut in interaction_sets[x]:
							interaction_sets[x][mut].append(decoy_energy)
						else:
							interaction_sets[x][mut] = [decoy_energy]

		for i, res in enumerate(interaction_sets):
			for v in res.values():
				v.sort(key=lambda x: (x[0][i], x[1]))

		return interaction_sets

	def pick_representative_deocys(self):
		"""
		Using favorable_interactions_by_pair_energy to create a starting list
		of protease mutations that yield favorable interactions with each 
		peptide residue, this function rearranges that list to organize the 
		decoys with favorable interactions by both peptide residue position, 
		and by what amino acid the peptide has in that position. It then trims
		the list of decoys representing each type of interaction to the three 
		with the most favorable interaction energy.
		"""
		interaction_sets = self.favorable_interactions_by_pair_energy()

		# Rearranging the interaction_sets
		rearranged_int_set = [{} for i in self.pep_res]
		for n, res in enumerate(interaction_sets):
			for k, v in res.items():
				for decoy in v:	
					aa = decoy[0][n]
					if aa not in rearranged_int_set[n]: 
					# First instance of this amino acid in this peptide space
						rearranged_int_set[n][aa]={k:[decoy[1:]]}
					else:
						if k not in rearranged_int_set[n][aa]:
						# First pairing of this mutation with this AA 
							rearranged_int_set[n][aa][k] = [decoy[1:]]
						else:
						# appending to existing mutation/aa pair
							rearranged_int_set[n][aa][k].append(decoy[1:])

		# Sorting the rearranged list
		for i in rearranged_int_set:
			for val in i.values():
				for v in val.values():
					v.sort(key=lambda x: x[0])
					while len(v) > 3: # Paring down to best three decoys
						v.pop(-1)

		return rearranged_int_set


def aggregated_report(name, mc_set):
	"""
	Whereas the mutations_by_seq_report lists the mutations each sequence 
	experienced, this report summarizes mutation occurrences across the entire
	available dataset.
	"""
	# Get full list of residues that could have been mutated
	ma = mutations_aggregate(mc_set)

	line_length = ma.output_length
	template = '{:12s}' * line_length + '\n'
	header = ['location', 'interacts', 'start_AA', 'start_type', 'mut_prob', 
				'final_AA', 'end_type', 'rel_prob', 'cleaved', 'pep_seq',  
				'seq_class', 'res_energy', 'int_energy', 'mut_prob_by_seq']

	with open(name, 'w') as r:
		r.write(template.format(*header))

		for line in ma.report:
			while len(line) < line_length:
				line.append('')
			r.write(template.format(*line))

	print name
	return ma


def representative_decoys_report(name, aggregate):
	"""
	Writes a file using the list of dictionaries generated by the 
	pick_representative_deocys in the mutations_aggregate object,
	tabulating the peptide position, the residue in that position, the 
	mutation producing a favorable interaction with that position, the best 
	energy from that interaction, and a command for PyMOL. The command will 
	delete everything currently in the viewer, then load the appropriate set 
	of PDB's exhibiting the interaction, run the compare_protease script, and
	zoom to the mutated residue in the protease. 
	"""
	position_ref = {0:'p6', 1:'p5', 2:'p4', 3:'p3', 4:'p2', 5:'p1'}

	rep_set = aggregate.representative_deocys
	header = ['Locus', 'Res', 'Mutant', 'Min_E', 'E_Change', 'Decoys']
	template = '{:12s}' * 5 + '{}\n'

	lines = []
	for n, i in enumerate(rep_set):
		locus = position_ref[n]
		for aa, v in i.items():
			for mut, decs in v.items():
				min_e = decs[0][0] # decoys are sorted, energy is first
				e_change = decs[0][1] # decoys are sorted, energy is first
				decoys = [x[2] for x in decs]
				decoys += [x[2].replace('designed', 'relaxed') for x in decs]
				res = str(mut.split('_')[0]) # Mutated residue on protease
				cmd = 'delete all; for pdb in ' + str(decoys) 
				cmd += ': cmd.load(pdb); @compare_protease.pml; '
				cmd += 'zoom res ' + res # Jump to relevant residue
				out_line = [locus, aa, mut, str(min_e), str(e_change), cmd]
				lines.append(out_line)

	lines.sort(key=lambda x: (x[0], x[2], x[1]))

	with open(name, 'w') as r:
		r.write(template.format(*header))
		for line in lines:
			r.write(template.format(*line))

	print name


##############################################################################
def main():
	# Getting user inputs
	args = parse_args()
	folder = args.folder
	base_name = basename(folder.rstrip('/'))

	# Intitialize Rosetta
	ros_opts = '-mute all'
	init(ros_opts)

	# Checking if analysis has already been run and pickled
	pickle_name = join(folder, base_name + '_mutations.pkl')

	if args.force_evaluate or not isfile(pickle_name):	
		# Making mutation_collection objects for each sequence
		report_files = sorted(glob(join(folder, '*.fasc')))
		for n, i in enumerate(report_files):
			if 'combined_reports.fasc' in i: # Skip file made by condense_fasc
				report_files.pop(n)

		# Simultaneously analyzing and writing to by-sequence report
		peptide_sequences = \
			[isolate_sequence_from_fasc(i) for i in report_files]
		seq_mutants = []
		report_name = join(folder, base_name + '_mutation_summary.txt')
		for n, i in enumerate(peptide_sequences):
			seq_muts = mutation_collection(folder, i, args.design_peptide)
			seq_mutants.append(seq_muts)

			if i == peptide_sequences[0]:
				mutations_by_seq_report(report_name, seq_muts, head=True)
			else:
				mutations_by_seq_report(report_name, seq_muts)

		# Pickling 
		with open(pickle_name, 'wb') as out:
			pickle.dump(seq_mutants, out, protocol=-1)

		print report_name

	else:
		# Retrieving pickled info
		with open(pickle_name, 'rb') as i:
			seq_mutants = pickle.load(i)

	# Making cross-sequence aggregated report
	agg_name = join(folder, base_name + '_by_mutations.txt')
	mutations_aggregate = aggregated_report(agg_name, seq_mutants)

	# Making list of representative decoys
	d_list_name = join(folder, base_name + '_look_at_me.txt')
	representative_decoys_report(d_list_name, mutations_aggregate)

if __name__ == '__main__':
	main()