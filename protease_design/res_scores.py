#!/usr/bin/python
from pyrosetta import *
from design_protease import res_to_design, apply_constraints
from glob import glob
from os.path import basename, join
opts='-mute core -mute protocols -mute basic -enzdes::cstfile ly104.cst -cst_fa_weight 1.0 -run:preserve_header'
init(opts)

cleaved = glob(join(basename('a_to_s_ly104_WT_CLEAVED_decoys'),"*best.pdb.gz"))
uncleaved = glob(join(basename('a_to_s_ly104_WT_UNcleaved_decoys'),"*best.pdb.gz"))

sf=create_score_function('ref2015_cst')
des_res=res_to_design('ly104_WT.pdb',8,[72,154])[2]
pep_res=range(197,208)

with open('res_enetgies.txt', 'w') as w:
	w.write('Cleaved\n')
	for i in cleaved:
		print i
		try:
			pose=apply_constraints(pose_from_pdb(i))
			sf(pose)
			energies=str(pose.energies()).split('\n')
			if i == cleaved[0]:
				w.write(energies[0].lstrip() + '\n')
			w.write(basename(i) + '\nNeighbor residues\n')
			for j in des_res:
				w.write(energies[j].lstrip() + '\n')
			w.write('Peptide residues\n')
			for j in pep_res:
				w.write(energies[j].lstrip() + '\n')
		except:
			continue

	w.write('\nUncleaved\n')
	for i in uncleaved:
		print i
		pose=apply_constraints(pose_from_pdb(i))
		sf(pose)
		energies=str(pose.energies()).split('\n')
		w.write(basename(i) + '\nNeighbor residues\n')
		for j in des_res:
			w.write(energies[j].lstrip() + '\n')
		w.write('Peptide residues\n')
		for j in pep_res:
			w.write(energies[j].lstrip() + '\n')
