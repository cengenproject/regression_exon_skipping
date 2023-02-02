
# Perform regression with different methods. Combination of explore_lasso_on_PSI_vs_TPM and explore_adaptive_lasso_on_PSI_vs_TPM

# Same as others, but for LASSO and adaptive LASSO, predict the deltaPSI instead of the PSI.
# What's more, try both with deltaPSI on natural scale, or taking the logit
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9


# copied from regressions_PSI_vs_TPM_deltaPSI.R to run this part on Ruddle with many replicates



# Inits ----

library(tidyverse)
library(glmnet)
library(furrr)

plan(multicore, workers = 10)


source("R/regression_functions.R")


out_file <- "data/intermediates/230202_regression_permutations_.qs"


# Read data ----


quantifs <- qs::qread("data/intermediates/230117_quantifs_filtered.qs")

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- qs::qread("data/intermediates/230117_tx_filtered.qs") |>
  filter(gene_id %in% putative_splice_factors)



events_coordinates <- read_tsv("data/export_for_arman/221111_events_coordinates.tsv")


# Filter events ----

# Keep only if enough reads supporting measure

quantifs_filtered_n_reads <- quantifs |>
  filter(nb_reads > 20)


# filter based on nb of samples that event was measured in

events_to_keep_n_samples <- quantifs_filtered_n_reads |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  filter(nb_samples > 120,
         nb_neurons > 32) |>
  pull(event_id)


quantifs_filtered_nsamples <- quantifs_filtered_n_reads |>
  filter(event_id %in% events_to_keep_n_samples)



# Filter events to remove those that are not DS between neuron types

events_to_keep_variability <- quantifs_filtered_nsamples |>
  filter(!is.na(PSI)) |>
  group_by(neuron_id, event_id) |>
  summarize(mean_PSI  = mean(PSI),
            sd_PSI = sd(PSI)) |>
  group_by(event_id) |>
  summarize(mean_PSI_btw_neurs  = mean(mean_PSI),
            sd_PSI_btw_neurs = sd(mean_PSI),
            .groups = 'drop') |>
  filter(sd_PSI_btw_neurs > 0.05) |>
  pull(event_id)




events_to_keep <- intersect(events_to_keep_n_samples,
                            events_to_keep_variability)

quantifs_filtered <- quantifs_filtered_n_reads |>
  filter(event_id %in% events_to_keep) |>
  filter(! is.na(PSI))

all.equal(sort(events_to_keep), sort(unique(quantifs_filtered$event_id)))



# Make SF expression as a matrix for use in regression

mat_sf_expression <- sf_expression |>
  select(transcript_id, sample_id, TPM) |>
  mutate(TPM = log(TPM + 1)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "TPM") |>
  column_to_rownames("sample_id") |>
  as.matrix() |>
  scale()
mat_sf_expression <- mat_sf_expression[,! apply(mat_sf_expression, 2, \(col) any(is.na(col)))]






####  Using permutation testing ----


regression_permutations <- expand_grid(event_id = sample(unique(quantifs_filtered$event_id), 3),
                                       method = c("lasso"),
                                       column = c("PSI"),
                                       shuffle = c(FALSE, rep(TRUE, 3))) |>
  mutate(res = future_pmap(list(event_id, method, column, shuffle),
                    \(event_id, method, column, shuffle) sparse_regression(event_id, method, column,
                                                                           shuffle, mat_sf_expression, quantifs_filtered),
                    .progress = TRUE),
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)

# qs::qsave(regression_permutations, out_file)


cat("Done, saved in ", out_file)
