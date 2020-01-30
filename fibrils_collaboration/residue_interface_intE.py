import argparse
from glob import glob
from os.path import basename, join
from pyrosetta import *
from pyrosetta.rosetta.core.simple_metrics.metrics import \
	TotalEnergyMetric, InteractionEnergyMetric
from pyrosetta.rosetta.core.select.residue_selector import \
	ChainSelector, ResidueIndexSelector

def parse_args():
	info = "Design a protease around a peptide sequence"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-d", "--directory", required=True,
		help="Pick a folder to analyze")
	args = parser.parse_args()
	return args

def calc_residue_interaction(pose, sf, res, int_target):
	ris = ResidueIndexSelector(str(res))
	inte = InteractionEnergyMetric()
	inte.set_scorefunction(sf)
	inte.set_residue_selectors(ris, int_target)
	return inte.calculate(pose)

args = parse_args()

init()
pdbs = glob(join(args.directory, '*.pdb'))
pdbs += glob(join(args.directory, '*.pdb.gz'))
pdbs.sort()
sf = get_fa_scorefxn()
tem = TotalEnergyMetric()
tem.set_scorefunction(sf)

csa = ChainSelector('A')
csb = ChainSelector('B')

report_name = basename(args.directory.rstrip('/')) + '_res_interaction_energies.csv'
print(args.directory)
print(report_name)
with open(join(args.directory, report_name), 'w') as w:
	for p in pdbs:
		print(p)
		pp = pose_from_pdb(p)
		sf(pp)
		chain_1_size = pp.split_by_chain()[1].total_residue()
		chain_2_size = pp.split_by_chain()[2].total_residue()
		this_res = [p, tem.calculate(pp)]
		for i in range(1, chain_1_size + 1):
			print(i)
			e = calc_residue_interaction(pp, sf, i, csb)
			this_res.append(e)
		print('chain change')
		for i in range(chain_1_size + 1, chain_1_size + chain_2_size + 1):
			print(i)
			e = calc_residue_interaction(pp, sf, i, csa)
			this_res.append(e)
		w.write(','.join([str(x) for x in this_res])+'\n')
