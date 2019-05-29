#!/bin/bash
for f in 5amu 5cpu 5dzu 5eu c35amu cyclo m5c; do com="python early_design_script.py ${f} \${job}"; python ../slurmit_BAY.py --job "${f}" --partition main --tasks 1 --requeue True --usearray True --array 1-25 --time 2:00:00 --begin now --outfiles fixed_designs/logs/${f}_job%a --command "$com"; done 