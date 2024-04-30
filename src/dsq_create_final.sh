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

dsq --job-file joblists/240415_final.txt --mem 15G --time 15:30:00 --partition day

