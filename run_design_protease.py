"""
Program for running presets on design_protease.py
To be populated with frequently-used proteases for convenience

Presently configured to be run from main khare_lab directory

Currently includes:
HCV
"""
import argparse

def parse_args():
	info = "Hand presets to design_protease.py, including which protease, \
	and which mutations"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-od", "--out_dir", required=True,
		help="Name an output directory for decoys")
	parser.add_argument("-name", "--name", type=str,
		help="How would you like to name your outputs? \
		(Default will use the name of the input PDB file.)")
	parser.add_argument("-seq", "--sequence", type=str,
		help="What substrate sequence do you want to thread")
	

def default_start_struct(protease):
	""" Pick out the starting PDB model for a given protease """
	if protease == 'HCV': return 'protease_design/start_proteases/HCV.pdb'

def default_cat_res(protease):
	""" Pick out the catalytic residues for a given protease """
	if protease == 'HCV': return '72 96 154'

def default_sequence(protease):
	""" Pick out default substrate sequence for a protease """
	if protease == 'HCV': return 'ADEMEECASHL'

def default_thread_site(protease, sequence):
	""" 
	For a given protease and sequence, pick out the appropriate threading site
	"""
	



