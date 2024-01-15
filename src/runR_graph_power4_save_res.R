#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=gr_pow4
#SBATCH -c 1
#SBATCH --mem-per-cpu=260G
#SBATCH --time=0:05:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


##do not use --mem=70G

echo "---------------  Running R script graph_power4_oncluster.R ---------------------"
echo "save results in more compact csv"
module load R
R --slave -f R/graph_power4_save_results.R
#R --slave -f R/graph_power_nosep_4.R
