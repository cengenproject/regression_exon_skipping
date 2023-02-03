#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=perm_test
#SBATCH -c 10
#SBATCH --mem=75G
#SBATCH --time=1-1:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "---------------  Running R script splice_graph.R ---------------------"

module load R
R --slave -f R/ruddle_permutation_test_PSI_vs_TPM.R

