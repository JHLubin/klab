#!/usr/bin/python
"""
Selects best scoring structures for design_protease simulations
"""
import argparse
from glob import glob
from os import mkdir
from os.path import basename, isdir, join
from shutil import copyfile
from sys import exit

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('source', type=str, help='Copy files from what directory?')
	parser.add_argument('destination', type=str, help='Copy files into what directory?')
	parser.add_argument('-n', '--n_models', default=10, type=int, 
		help='Pick how many decoys? (Default: 10)')
	parser.add_argument('-t', '--score_type', default='total_score', type=str, 
		help='Select models based on which score? (Default: total_score)')
	args = parser.parse_args()
	return args


def collect_fasc_files(directory):
	""" 
	Collects all .fasc files in a directory, excluding combined files made by
	condense_fasc.py.
	"""
	folder_search = join(directory, "*.fasc")
	fasc_files = glob(folder_search)
	fasc_files.sort()

	# Prevent self-reference if re-running
	for n, i in enumerate(fasc_files):
		if 'combined_reports' in i:
			fasc_files.pop(n)

	return fasc_files


def sort_score_file(score_file, score_term, n_models):
	""" Reads a score file and sorts by a given column """
	with open(score_file, 'r') as r:
		readout = r.readlines()[1:]

	scores_dict = {}
	for line in readout:
		linesplit = line.replace(':','').split()
		decoy = None
		score = None

		# Getting file name and selected score
		for i, j in enumerate(linesplit):
			if j == 'filename':
				decoy = basename(linesplit[i + 1])
			if j == score_term:
				score = float(linesplit[i + 1])

		# Throwing errors of file name or score are not found
		if decoy == None:
			print "Not finding decoy name column correctly; check score file."
			exit(1)

		if score == None:
			print "Could not find specified score term; check input/score file."
			exit(1)

		# Checking for duplicate names, and if none, adding to dict
		if decoy in scores_dict:
			print decoy, 'appeared a second time. Check source files.'
			exit(1)
		else:
			scores_dict[decoy] = [score, line]

	# Sorting, filtering
	sorted_scores_dict = sorted(scores_dict.items(), key=lambda x: x[1][0])
	selection = sorted_scores_dict[:n_models]
	sorted_files = [i[0] for i in selection]
	sorted_lines = [i[1][1] for i in selection]

	return sorted_files, sorted_lines


def main():
	# Getting user inputs
	args = parse_args()

	# Getting fasc files from source directory
	source = args.source
	fasc_files = collect_fasc_files(source)

	# Establishing output directory
	dest = args.destination
	if not isdir(dest):
		mkdir(dest)

	files_to_copy = []
	for sf in fasc_files:
		# Sorting score file
		decoys, scorelines = sort_score_file(sf, args.score_type, args.n_models)

		# Adding PDB files to list for to be copied over
		files_to_copy += decoys

		# Creating shortened score file
		new_sf = join(dest, basename(sf))
		with open(new_sf, 'w') as w:
			w.write('\n') # Takes place of header line
			for l in scorelines:
				w.write(l)

		# Copying over threaded model
		threaded_pdb = basename(sf).replace('_designed.fasc', '.pdb.gz')
		copyfile(join(source, threaded_pdb), join(dest, threaded_pdb))

	# Copying over decoys
	print "Copying files:"
	for fi in files_to_copy:
		print '\t', fi
		designed = fi
		relaxed = fi.replace('designed', 'relaxed')
		copyfile(join(source, relaxed), join(dest, relaxed))
		copyfile(join(source, designed), join(dest, designed))

if __name__ == '__main__':
	main()