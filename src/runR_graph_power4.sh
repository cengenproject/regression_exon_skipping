#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=gr_pow4
#SBATCH -c 20
#SBATCH --mem=2G
#SBATCH --time=3:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "---------------  Running R script graph_power4_oncluster.R ---------------------"

module load R
R --slave -f R/graph_power4_oncluster.R
