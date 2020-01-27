from glob import glob
from pyrosetta import *
from pyrosetta.rosetta.core.simple_metrics.metrics import TotalEnergyMetric, InteractionEnergyMetric
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import PerResidueEnergyMetric
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, ResidueIndexSelector
init()
pdbs = glob('*.pdb')
pose=pose_from_pdb('htra1_pdz_asyn_ext_1_19.pdb')
sf=get_fa_scorefxn()
prem = PerResidueEnergyMetric()
prem.set_scorefunction(sf)
tem = TotalEnergyMetric()
tem.set_scorefunction(sf)

results = []
for p in pdbs:
    pp = pose_from_pdb(p)
    results.append([p, tem.calculate(pp), prem.calculate(pp)])
s=prem.calculate(pose)

with open('res_results.txt','w') as w:
    for r in results:
        res_scores = []
        for i in range(1,116):
            res_scores.append(r[2][i])
        ol = [r[0], r[1]] + res_scores
        w.write(','.join([str(x) for x in ol])+'\n')

csb = ChainSelector('B')
with open('res_interactions.csv','w') as w:
	for p in pdbs:
	    pp = pose_from_pdb(p)
	    tot_e = tem.calculate(pp)
	    this_res = [p, tot_e]
	    for i in range(1,115):
	    	ris = ResidueIndexSelector(str(i))
	    	inte = InteractionEnergyMetric()
	    	inte.set_scorefunction(sf)
	    	inte.set_residue_selectors(ris, csb)
	    	this_res.append(inte.calculate(pp))
	    w.write(','.join([str(x) for x in this_res])+'\n')