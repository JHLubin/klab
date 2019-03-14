#!/usr/bin/python
"""
Combines FASC files from PyRosetta job distributor and converts them into a 
tabular form. 
"""
import argparse
from glob import glob
from os.path import basename, join

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("directory", type=str, 
		help="Read fasc files from what directory?")
	parser.add_argument("-z", "--zipped", action="store_true", default=False,
		help="Add .gz suffix to decoys?")
	parser.add_argument("-kf", "--keep_folder", action="store_true", 
		default=False, help="Keep the folder name of the decoy? (By default, \
		only the basename will appear.)")
	args = parser.parse_args()
	return args


def text_to_numbers(in_list):
	""" 
	Takes a list that includes numbers and text all stored as strings.
	Returns the same list with the numbers converted to floats.
	"""
	new_list = []
	for i in in_list:
		try:
			new_list.append(float(i))
		except ValueError:
			new_list.append(i)

	return new_list


def cleanup_mutations_section(report_lines, start_point):
	""" 
	Some sets are missing mutable residues present in others. This function 
	maintiains column alignment across all sets by inserting what will appear
	as blank cells where a mutable residue is missing.
	"""
	max_len = max([len(line) for line in report_lines])
	res_columns = list(range(start_point, max_len, 4))

	for c in res_columns:
		mutated_residues = [int(line[c]) for line in report_lines if line[c] != "NONE"]
		if len(mutated_residues) == 0:
			return
		first_des_res = min(mutated_residues)
		for line in report_lines:
			if int(line[c]) != first_des_res:
				for i in range(4):
					line.insert(c, '=""')

	return


def main():
	# Getting user inputs
	args = parse_args()

	# Getting fasc files
	folder = args.directory
	base_name = basename(folder.rstrip('/'))
	folder_search = join(folder, "*.fasc")
	fasc_files = glob(folder_search)
	fasc_files.sort()
	# Prevent self-reference if re-running
	for n, i in enumerate(fasc_files):
		if 'combined_reports' in i:
			fasc_files.pop(n)

	# Collecting fasc lines
	headline = []
	report_lines = []
	mut_section_ind = 0
	for f in fasc_files:
		# Reading in individual FASC file
		with open(f, 'r') as read:
			f_lines = read.readlines()
			f_lines.pop(0) # First line is not useful
			lines_data = []
			for i in f_lines:
				if ' mutations:' in i: # Specific to design_protease
					start_mutations = i.find(' mutations:') 
					scores_section = i[:start_mutations + 1].split()
					prot_mutations_section = i[start_mutations + 1:].split()
					line_data = scores_section[1::2] + prot_mutations_section
				else:
					scores_section = i.split()
					line_data = scores_section[1::2]

				if i == f_lines[0] and f == fasc_files[0]:
					headline = scores_section[::2]
					mut_section_ind = int(len(scores_section) / 2 + 1) # Specific to design_protease

				lines_data.append(text_to_numbers(line_data))

			lines_data.sort(key=lambda x: x[1]) # Sorting by total score
			report_lines += lines_data

	# Making combined report
	report_name = join(folder, base_name + '_combined_reports.fasc')

	with open(report_name, 'w') as r:
		# Making template and header
		head_length = len(headline)
		template = '{:50}' + '{:25}' * (head_length - 1)
		r.write(template.format(*headline) + '\n')

		# Adding in lines
		cleanup_mutations_section(report_lines, mut_section_ind)
		for line in report_lines:
			if not args.keep_folder: # Stripping decoy name to base
				line[0] = basename(line[0])
			if args.zipped: # Adding zip suffix if that option is used
				line[0] += '.gz'
			line_out = template.format(*line[:head_length])
			line_out += '   '.join(line[head_length + 1:])
			r.write(line_out + '\n')

		print(report_name)

if __name__ == '__main__':
	main()