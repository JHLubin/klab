#!/usr/bin/python
""" 
PyRosetta4, Python 2.7
Joseph Lubin, 2018
"""
from pyrosetta import *
from glob import glob
from os.path import basename, join
init('-mute all')

def get_original_seq(pdb_name):
	""" Gets 6-letter peptide sequence from name of PDB file """
	split = pdb_name.split('_')
	long_seq = split[3]
	short_seq = long_seq[1:7]

	return short_seq


def is_cleaved(seq, cleaveds, uncleaveds):
	""" Checks whether a sequence is on the cleaved or uncleaved list """
	if seq in cleaveds:
		return 'cleaved'
	elif seq in uncleaveds:
		return 'uncleaved'
	else: 
		return 'new'


# Getting PDB files
folder = 'new_design_pep+prot_fastdesign'
designed_search = '*designed*.pdb'
pdb_set = sorted(glob(join(folder, designed_search)))

# Getting lists of cleaved/uncleaved sequences
file_cleaved = 'list_cleaved.txt'
with open(file_cleaved, 'r') as fc:
	cread = fc.readlines()
	list_cleaved = [i.strip() for i in cread]

file_uncleaved = 'list_uncleaved.txt'
with open(file_uncleaved, 'r') as fu:
	uread = fu.readlines()
	list_uncleaved = [i.strip() for i in uread]

# Starting report
report = 'new_design_pep+prot_fastdesign/mutated_peptides.txt'
template = '{:50s}' + '{:12s}' * 4 + '{:20s}\n'
header = ['PDB', 'NewSeq', 'Cleaved?', 'OldSeq', 'Cleaved?', 'MutableResidues']
with open(report, 'w') as out:
	out.write(template.format(*header)) 


# Reading PDB file for sequence
aggregated_list = {}
for pdb in pdb_set:
	pbase = basename(pdb)
	original_seq = get_original_seq(pbase)
	original_status = is_cleaved(original_seq, list_cleaved, list_uncleaved)
	pose = pose_from_pdb(pdb)
	pose_sequence = pose.sequence()
	designed_seq = pose_sequence[197:203]
	designed_status = is_cleaved(original_seq, list_cleaved, list_uncleaved)

	des_res = [58, 70, 138, 147, 150, 152, 153, 169, 170, 171, 172, 173, 174, 175, 176, 177, 183]
	res_list = ''
	for r in des_res:
		res_list += pose_sequence[r-1]

	# Writing to output file
	outline = [pbase, designed_seq, designed_status, original_seq, original_status, res_list]
	with open(report, 'a') as out:
		line = template.format(*outline)
		print pbase, line[:-2]
		out.write(line)

	# Collecting by mutated sequence
	if designed_seq in aggregated_list:
		aggregated_list[designed_seq].append(pbase)
	else:
		aggregated_list[designed_seq] = [pbase]

# Writing output by designed sequence
group_summary = 'new_design_pep+prot_fastdesign/result_pep_sequences.txt'
with open(group_summary, 'w') as go:
for s, pl in aggregated_list.iteritems():
go.write('\nSequence:\t' + s + '\n')
for p in pl:
go.write('\t\t' + p + '\n')

"""
with open('new_design_pep+prot_fastdesign/mutated_peptides.txt','r') as r:
    lines = r.readlines()

aggregated_list = {}

for line in lines:
    s=line.split()
    seq=s[1]
    if seq in aggregated_list:
        aggregated_list[seq].append(s[0])
    else:
        aggregated_list[seq]=[s[0]]

lengths = {}

for n, i in aggregated_list.iteritems():
    if len(i) in lengths:
        lengths[len(i)]+=1
    else:
        lengths[len(i)]=1

lengths = {1: 121, 2: 39, 3: 22, 4: 9, 5: 3, 6: 4, 7: 6, 8: 2, 9: 2, 
10: 4, 12: 1, 15: 1, 17: 1, 21: 1}
"""