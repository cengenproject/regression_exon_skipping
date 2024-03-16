# Regressions to find regulators of splicing



# Graph power 4

Main file: `graph_power4_oncluster.R`. Can be called from CLI with right arguments: use `src/dsq_create.sh` to create a job array. This file describes a workflow using the steps defined in `functions_steps.R`, and there are many helper functions in `loss_functions.R` (which could be a package).


Note: for permutations, if `permutation` == 0 , no permutation, any number bigger than 0 the value doesn't matter. The seed for the entire script is set by `sum`(permutations)` (to avoid having independent runs identical), but there shouldn't be much else affected (the diagonal in CLIME for a big penalty is also randomly generated).


After being run, download all csv files, and use `analysis_from_cluster.R`, which sources `analysis_helpers.R`.



# Older versions

## "Forward" fit

In `regression_functions.R`, fit PSI = C * TPM

## Simulations

In `simulation_perm_test.R`

## Permutation tests

Adapted to run on Ruddle: `ruddle_permutation_test_PSI_vs_TPM.R`, `ruddle_permutation_test_PSI_vs_TPM_sim.R`, then results interpretation in `simulation_perm_test.R` and `results_permutation_test.R`.


