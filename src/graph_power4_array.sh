#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=gr_pow4
#SBATCH -c 1
#SBATCH --mem=25G
#SBATCH --time=13:40:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


##do not use --mem=70G

echo "---------------  Running R script graph_power4_oncluster.R ---------------------"
echo "240126_params_noperm_7penalties"
module load R
Rscript R/graph_power4_oncluster.R \
  --date '240315' \
  --exonsInput 'PSI' \
  --transformation 'npnshrink' \
  --permutations '0:20' \
  --penalties 'c(10,5,1,.7,.6,.5,.4,.3,.2,.1,.05)' \
  --algo 'QUIC'

#R --slave -f R/graph_power_nosep_4.R
