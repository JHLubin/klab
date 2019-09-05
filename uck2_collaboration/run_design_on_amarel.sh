#!/bin/bash
for f in 5amu 5cpu 5dzu 5eu c35amu cyclo m5c
 do com="python early_design_script.py ${f} \${job}"
 python ../slurmit_BAY.py --job "${f}" --partition main --tasks 1 --requeue True --usearray True --array 1-25 --time 2:00:00 --begin now --outfiles fixed_designs/logs/${f}_job%a --command "$com"
 done 


odir=/scratch/jhl133/uck2_collaboration/designing_enzyme/all_ala; mkdir -p $odir/logs; for ba in c u; do for va in 5am 5cp 5dz 5e c35am cyclo m5 wt; do jobname=$va$ba; com="python uck2_collaboration/early_design_script.py -od $odir -l $va$ba -n 10 -x \${jobname}_\${job} -ala"; python slurmit_BAY.py --job "${jobname}" --partition main --tasks 1 --usearray True --array 1-10 --requeue True --time 3-00:00:00 --begin now --outfiles $odir/logs/${job}_job%a --command "$com"; done; done
odir=/scratch/jhl133/uck2_collaboration/designing_enzyme/wt_start; mkdir -p $odir/logs; for ba in c u; do for va in 5am 5cp 5dz 5e c35am cyclo m5 wt; do jobname=$va$baw; com="python uck2_collaboration/early_design_script.py -od $odir -l $va$ba -n 10 -x \${job}"; python slurmit_BAY.py --job "${jobname}" --partition main --tasks 1 --usearray True --array 1-10 --requeue True --time 3-00:00:00 --begin now --outfiles $odir/logs/${job}_job%a --command "$com"; done; done

odir=/scratch/jhl133/uck2_collaboration/designing_enzyme_round_2; mkdir -p $odir/logs; for ba in c u; do for va in 5am 5cp 5dz 5e c35am cyclo m5 wt; do jobname=a$va$ba; com="python uck2_collaboration/early_design_script.py -od $odir -l $va$ba -n 5 -x all_ala_\${jobname}_\${job} -ala"; echo python slurmit_BAY.py --job "${jobname}" --partition main --tasks 1 --usearray True --array 1-5 --requeue True --time 3-00:00:00 --begin now --outfiles $odir/logs/${job}_job%a --command "$com"; done; done; for ba in c u; do for va in 5am 5cp 5dz 5e c35am cyclo m5 wt; do jobname=w$va$baw; com="python uck2_collaboration/early_design_script.py -od $odir -l $va$ba -n 5 -x wt_start\${job}"; echo python slurmit_BAY.py --job "${jobname}" --partition main --tasks 1 --usearray True --array 1-5 --requeue True --time 3-00:00:00 --begin now --outfiles $odir/logs/${job}_job%a --command "$com"; done; done
