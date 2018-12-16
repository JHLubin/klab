#!/bin/sh
#SBATCH --job-name=manaseq
#SBATCH --partition=p_sdk94_1,main
#SBATCH --requeue
#SBATCH --ntasks-per-node=28
#SBATCH --mem=128GB
#SBATCH --nodes=2
#SBATCH --time=1:0:0
#SBATCH --array=0-55
#SBATCH --output manasi_sequences/outerr/esd.%j.%N.out
#SBATCH --error manasi_sequences/outerr/esd.%j.%N.err
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhl133@scarletmail.rutgers.edu

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

python design_protease.py -s a_to_s_ly104_WT.pdb -od manasi_sequences -cseq pep_seq_lists/manasi_cleaved.txt -useq pep_seq_lists/manasi_uncleaved.txt -cons ly104.cst -d 20 > manasi_sequences/outerr/manasi_sequences_${SLURM_ARRAY_TASK_ID}.out
