# Regressions to find regulators of splicing



# Graph power 4


## Data preprocessing

The PSIs and TPMs are first evaluated elsewhere (repos `quantif_exon_skipping` and `stringtie_quantif`).

`R/prefilter_quantifs.R` takes PSI and TPMs and prefilters to keep only relevant cases. `R/prep_graph_power4.R` makes the train/test split, and saves the matrices ready to be used directly in main script.


## Network computation

Main file: `graph_power4_oncluster.R`. Can be called from CLI with right arguments: use `src/dsq_create.sh` to create a job array. This file describes a workflow using the steps defined in `functions_steps.R`, and there are many helper functions in `loss_functions.R` (which could be a package).


Note: for permutations, if `permutation` == 0 , no permutation, any number bigger than 0 indicates that permutation should be performed (the value doesn't matter). The seed for the entire script is set by `sum`(permutations)` to avoid having independent runs identical, but there shouldn't be much else affected (the diagonal in CLIME for a big penalty is also randomly generated).

### Runs

Everything run as a dSQ (dead simple queue) joblist. In `src/dsq_create_*`, scripts to create the job from a joblist in `joblists/`.

* First run joblist `noperm_glassoQuicClime.txt` to try all algorithms with a few penalties (used to generate heatmap Fig. 7B). In parallel, run `noperm_noperm_Scio.txt` with less memory, more time (note the names are wrong, in fact it's CLIME that is performed separately).

* Then run joblist `nopermglasso_k.txt`, focusing on a single algorithm but trying out different k in knn imputation (Fig. S5B).

* Then "perm.txt" uses single algorithm and k, and runs permutation tests at all penalties (Fig. S5C, Fig. 7C). Finally, "final.txt" reruns a single time at best penalty (and algo and k, and a single permutation), but saves the full intermediate object for in-depth examination.

## Post-analysis

After being run, download all csv files, and use `analysis_from_cluster.R`, which sources `analysis_helpers.R` and performs the comparisons (Fig 7B,C, S5B,C). The adjacency matrix is loaded and analyzed in `analysis_network.R` to create the subnetworks (Fig 7D,E, S5D).

Some "manual" explorations are in `graph_power4_postanalysis.R` (used to e.g. compare the impact of individual changes, not included in the paper).



# Older versions

## "Forward" fit

In `regression_functions.R`, fit PSI = C * TPM

## Simulations

In `simulation_perm_test.R`

## Permutation tests

Adapted to run on Ruddle: `ruddle_permutation_test_PSI_vs_TPM.R`, `ruddle_permutation_test_PSI_vs_TPM_sim.R`, then results interpretation in `simulation_perm_test.R` and `results_permutation_test.R`.


