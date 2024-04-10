#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=dsq_create
#SBATCH -c 1
#SBATCH --mem=100M
#SBATCH --time=00:02:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

ml dSQ

dsq --job-file src/240410_noperm_glasso_k.txt --mem 20G --time 2:30:00

