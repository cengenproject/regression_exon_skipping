#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=gr_pow4
#SBATCH -c 6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=5-3:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


##do not use --mem=70G

echo "---------------  Running R script graph_power4_oncluster.R ---------------------"
echo "50 permutations, separate Nincl/excl"
module load R
R --slave -f R/graph_power4_oncluster.R
#R --slave -f R/graph_power_nosep_4.R