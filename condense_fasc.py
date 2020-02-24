#!/usr/bin/python
import argparse
from glob import glob
import json
import os 
from shutil import copyfile

def parse_args():
	info = "Combines FASC files from PyRosetta job distributor and converts \
	them into a more convenient csv table. Also extracts best decoys by a \
	given criterion if desired."
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("directory", type=str, 
		help="Read .fasc files from what directory?")
	parser.add_argument("-n", "--num_decoys", type=int, default=0, 
		help="Collect the top [how many] decoys? By default, no top decoy \
		selection will be taken.")
	parser.add_argument("-od", "--out_directory", type=str,  
		help="Where should best decoys be copied? If num_decoys is given and \
		there is no input for this option, by default a subfolder will be \
		created in the given directory called \
		[original_dir]_best_[number]_decoys.")
	parser.add_argument("-s", "--selection_criterion", type=str, 
		default='total_score', help="By what property should best decoys be \
		selected? Default is total_score.")
	parser.add_argument("-m", "--merge_fascs", action='store_true', 
		help="If multiple .fasc files were created for a single set, this \
		option allows them to be sorted together. It is assumed that the \
		difference will only be in numbers after an underscore, and that the \
		beginning name will be the same for all clustered .fasc files.")
	args = parser.parse_args()
	return args


def fasc_to_dict_list(fasc_file, sortby='total_score'):
	"""
	Reads in a fasc file and converts it to a list of dicts. Then sorts the list 
	by a given criterion ('total_score' by default) and returns it.
	Note: does not work for files from older versions of Rosetta, which output 
	text strings that don't resemble dicts. Also adds several items to the decoy 
	dict: basename, set (which .fasc file or group the decoy was from), and rank
	(based on the sorting criterion). Rank is one-indexed, not zero-indexed. 
	"""
	fasc_lines = []
	# Reading in FASC files
	if isinstance(fasc_file, list):
		fasc_base = os.path.basename(fasc_file[0])
		fasc_cluster = fasc_base[:fasc_base.rfind('_')]
		for f in fasc_file:
			with open(f, 'r') as read:
				fasc_lines += read.readlines()
	else:
		fasc_base = os.path.basename(fasc_file)
		fasc_cluster = fasc_base[:fasc_base.rfind('.fasc')]
		with open(fasc_file, 'r') as read:
			fasc_lines += read.readlines()

	# Collecting decoy data from FASC file as a list of dicts
	decoys_data = [json.loads(fl) for fl in fasc_lines]

	# Sorting list
	decoys_data.sort(key=lambda x: x[sortby])

	# Add basename, set, and rank to each decoy
	for n, i in enumerate(decoys_data):
		dec_base = os.path.basename(i['filename'])
		i['decoy_basename'] = dec_base[:dec_base.index('.')] # Remove file extension
		i['decoy_set'] = fasc_cluster
		i['decoy_rank'] = n + 1 

	return decoys_data


def make_condensed_fasc_header(decoy_list):
	""" 
	For a given list of dicts representing decoys, creates a list for to be 
	used as a header. Collects all keys from the dicts in the decoy_list, 
	including any that are not present in all decoys. Eliminates pdb_name and 
	decoy from the list, since they are redundant, and nstruct since that is 
	unnecessary when the list is converted to CSV. Moves filename, decoy_set, 
	and decoy_rank to the start of the list, followed by total_score and 
	anything else with 'total', followed by anything with 'constraint'. 
	Preappends a column for the decoy basename followed by a variable number of 
	columns for non-identical underscore-separated parts of the decoy basename
	so it is easy to sort decoys later. Makes a masking list for the decoys' 
	names so that splitting them later, they can be appropriately tabulated with 
	the non-identical parts included and the identical parts excluded.
	"""
	# Initialize header
	header = ['decoy_basename']
	keylist = []

	# Collect all keys
	for decoy in decoy_list:
		for k in decoy.keys():
			if k not in keylist:
				keylist.append(k)
	keylist.remove('pdb_name')
	keylist.remove('decoy')
	keylist.remove('nstruct')
	keylist.sort()

	# Collect list of basenames split by underscore
	basenames = [d['decoy_basename'].split('_') for d in decoy_list]

	# Add columns to header for decoy name differences, make name mask
	decoy_name_mask = []
	for i in zip(*basenames):
		if len(set(i)) > 1: # True if not all elements are identical
			header.append('')
			decoy_name_mask.append(1)
		else:
			decoy_name_mask.append(0)

	# Moving properties from keylist to header
	header += ['filename', 'decoy_set', 'decoy_rank', 'total_score']

	# Moving any totals from keylist to header
	header += [k for k in keylist if 'total' in k and k not in header]

	# Moving any constraints from keylist to header
	header += [k for k in keylist if 'constraint' in k and k not in header]

	# Moving everything else from keylist to header
	header += [k for k in keylist if k not in header]

	return header, decoy_name_mask


def convert_dec_dict_to_csv_line(decoy_dict, header, mask):
	"""
	Takes in a decoy as a dict and populates a list of the values of that dict  
	in order matching the given header, filling in blanks if the decoy does not 
	have a given key. Uses the mask to populate columns with non-identical decoy 
	name features. Then converts the list to a comma-separated string.
	"""
	# Collect basename to initialize output line
	outlist = [decoy_dict['decoy_basename']]

	# Add non-identical name components to output line
	basename_split = decoy_dict['decoy_basename'].split('_')
	outlist += [i for i,j in zip(basename_split, mask) if j]

	# Add all other values in order according to header
	for score in header[header.index('filename'):]:
		if score in decoy_dict.keys():
			outlist.append(decoy_dict[score])
		else:
			outlist.append('')

	# Convert to string
	return ','.join([str(i) for i in outlist])


def check_make_folder(directory):
	"""
	Checks whether a directory exists, and creates it if it doesn't
	"""
	if not os.path.isdir(directory):
		os.makedirs(directory)

	return


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


def main(args):
	# Getting fasc files
	folder = args.directory.rstrip('/')
	base_name = os.path.basename(folder)
	folder_search = os.path.join(folder, "*.fasc")
	fasc_files = glob(folder_search)
	fasc_files.sort()

	# Group .fasc files if flag is given
	if args.merge_fascs:
		fasc_groups = {}
		for f in fasc_files:
			fasc_base = os.path.basename(f)
			fasc_cluster = fasc_base[:fasc_base.rfind('_')]
			if fasc_cluster in fasc_groups:
				fasc_groups[fasc_cluster].append(f)
			else: 
				fasc_groups[fasc_cluster] = [f]
		fasc_files = list(fasc_groups.values())

	# Collecting fasc lines
	master_decoy_list = []
	for f in fasc_files:
		fasc_extract = fasc_to_dict_list(f, sortby=args.selection_criterion)
		master_decoy_list += fasc_extract

	# Making header line and decoy name mask
	header, name_mask = make_condensed_fasc_header(master_decoy_list)

	# Making combined report
	report_name = os.path.join(folder, base_name + '_combined_reports.csv')
	with open(report_name, 'w') as r:
		# Writing header
		r.write(', '.join(header) + '\n')

		# Writing decoy lines
		for decoy in master_decoy_list:
			dec_line = convert_dec_dict_to_csv_line(decoy, header, name_mask)
			r.write(dec_line + '\n')

	print(report_name)

	# Copying best decoys
	if args.num_decoys:
		# Name output directory
		if args.out_directory:
			outdir = args.out_directory
		else:
			foldername = '{}_best_{}_decoys'.format(base_name, args.num_decoys)
			outdir = os.path.join(folder, foldername)
		check_make_folder(outdir)

		# Collect best decoys
		decoys_to_collect = []
		for decoy in master_decoy_list:
			if decoy['decoy_rank'] <= args.num_decoys:
				# Make name independent of file extension
				decoyname = os.path.join(folder, 
					decoy['decoy_basename'] + '.*')
				decoys_to_collect += glob(decoyname)

		# Copy best decoys
		for decoy in decoys_to_collect:
			print(os.path.basename(decoy))
			decoy_out = decoy.replace(folder, outdir)
			copyfile(decoy, decoy_out)


if __name__ == '__main__':
	args = parse_args()
	main(args)