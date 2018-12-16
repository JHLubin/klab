#!/usr/bin/python
"""
	No, this script won't return a bunch of models in that sense.
	However, this script can return to you the names of the top x models by score. Or even
	copy those models over to a directory and reorder them based on score. It all depends
	on the flags you give it. At a minimum, please specify the directory with the models
	and the score file.

	Example usage:
		python get_top_models.py ./pdbs/ ./score.sc --print_only
"""
import argparse
import os
import sys
import shutil

def parse_args():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('directory', type=is_dir, help='a directory of pdb files')
	parser.add_argument('score_file', type=argparse.FileType('r'), help='a score file corresponding to the pdbs')
	parser.add_argument('-n', '--n_models', default=10, type=int, help='the number of models to return, default is 10')
	parser.add_argument('-x', '--model_extension', default='.pdb', type=str, help='the extension of the models, default is .pdb')
	parser.add_argument('-c', '--column', default='total_score', type=str, help='the column to sort by, default is total_score')

	option_group = parser.add_mutually_exclusive_group(required=True)
	option_group.add_argument('--print_only', action='store_true', help='give this flag to print top n models')
	option_group.add_argument('--store', action='store', metavar='output_directory',type=is_dir, help='give this flag to copy over the top n models to the sepcified directory')
	args = parser.parse_args()
	return args

def is_dir(potential_dir):
	""" Quick and dirty function to check if input is a dir. """
	if not os.path.isdir(potential_dir):
		raise argparse.ArgumentTypeError("{0} is not a directory.".format(potential_dir))
	else:
		return potential_dir

def sort_score(score_file,column,n_models):
	""" give a score file (opened in r mode), parse the top models by overall score
		returns a list of sorted tuples [(model_name, score),...]
	""" 
	header = score_file.readline().split()
	score_column = None
	model_name_column = None

	for i in range(len(header)):
		if header[i] == column:
			score_column = i

		if header[i] == 'description':
			model_name_column = i

	if score_column is None:
		print 'Could not find the specified score column ({0}) in score file, sorry!'.format(column)
		sys.exit(1)

	if model_name_column is None:
		print 'Could not find description column in score file, sorry!'
		sys.exit(1)

	#initialize score_dict = {model:score} and store score file info
	score_dict = {}

	for line in score_file:
		try:
			score_dict[line.split()[model_name_column]] = float(line.split()[score_column])
		except KeyError:
			print 'Your score file has 2 models with the same name... why?'
			raise

	sorted_models = sorted(score_dict.items(), key=lambda x: x[1])

	return sorted_models[0:n_models]

def main():
	args = parse_args()

	sorted_models = sort_score(args.score_file,args.column,args.n_models)

	#getting proper tab formatting based on len of model name
	print 'model{0}score'.format('\t'*(len(sorted_models[0][0])/5-1))
	for model in sorted_models:
		print '{0}\t{1}'.format(model[0],model[1])

	if args.print_only == False:
		print 'Copying models to {0}'.format(args.store)
		print 'This will overwrite...'
		model_counter = 1
		for model in sorted_models:
			model_name = '{0}{1}'.format(model[0],args.model_extension)
			model_path = os.path.abspath(os.path.join(args.directory,model_name))

			modified_name = 'top{0}_{1}'.format(model_counter,model_name)
			model_counter += 1

			output_path = os.path.abspath(os.path.join(args.store,modified_name))
			shutil.copyfile(model_path,output_path)

if __name__ == "__main__":
	main()