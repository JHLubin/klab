#!/usr/bin/python
"""
When downloading a new PDB file, relax with coordinate constraints to eliminate clashes.

Requires a PDB file input.

Options:
Name (-n, string): change the output PDB name from [original_name]_relaxed.pdb
Score function (-sf, string): change the score function from the default of ref2015_cst
Catalytic residues (-cat, int, multiple accepted): list residues that should not be moved 

"""
from pyrosetta import *
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.relax import FastRelax
from os.path import basename

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("pdb_file", help="What PDB file do you want to relax?")
	parser.add_argument('-n', "--name", 
		help="What do you want to name the relaxed PDB? (Default appends '_relaxed' to original name.)")
	parser.add_argument('-sf', "--score_function", default='ref2015_cst',
		help="What score function do you want to use? (Default: ref2015_cst)")
	parser.add_argument('-cat', "--cat_res", type=int, nargs='+',
		help="List residues that should be immobile, separated by spaces")
	args = parser.parse_args()
	return args


def main(args):
	# Creating coordinate constraints for the entire molecule
	cg = CoordinateConstraintGenerator()
	ac = AddConstraints()
	ac.add_generator(cg)

	# Creating FastRelax protocol with the given score function
	fr = FastRelax()
	sf = create_score_function(args.score_function)
	fr.set_scorefxn(sf)

	# Creating a movemap with backbone fixed, side chains mobile
	mm = MoveMap()
	mm.set_bb(False)
	mm.set_chi(True)
	if args.cat_res:
		# Side chains fixed for specified catalytic residues
		for i in args.cat_res:
		    mm.set_chi(i, False)
	fr.set_movemap(mm)

	# Loading PDB file, applying constraints, relaxing
	pose = pose_from_pdb(args.pdb_file)
	ac.apply(pose)
	fr.apply(pose)




if __name__ == '__main__':
	args = parse_args()
	init('-cst_fa_weight 1.0')
	main(args)

# '-relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains'