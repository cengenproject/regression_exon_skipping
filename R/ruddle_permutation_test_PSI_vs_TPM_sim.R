
# Perform regression with different methods. Combination of explore_lasso_on_PSI_vs_TPM and explore_adaptive_lasso_on_PSI_vs_TPM

# Same as others, but for LASSO and adaptive LASSO, predict the deltaPSI instead of the PSI.
# What's more, try both with deltaPSI on natural scale, or taking the logit
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9


# copied from regressions_PSI_vs_TPM_deltaPSI.R to run this part on Ruddle with many replicates
# Run on results of simulation, ignore some columns that have become meaningless


cat("Starting at ",date(),"\n")

# Inits ----

library(tidyverse)
library(glmnet)
library(furrr)

plan(multicore, workers = 10)


source("R/regression_functions.R")


out_file <- "data/intermediates/230203_sim1_regression_permutations_psi.qs"


# Read data ----


quantifs <- qs::qread("data/intermediates/simultation/230206_quantifs_simulation1.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_sf_simulation1.qs")





# Events already filtered

quantifs_filtered <- quantifs



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


# Predict DeltaPSI instead of PSI
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9 for rationale

logit <- function(x){
  stopifnot(all(x >= 0 & x <= 1))
  x[x == 0] <- .01
  x[x == 1] <- .99
  log(x/(1-x))
}

quantifs_filtered <- quantifs_filtered |>
  group_by(event_id) |>
  mutate(dPSI_nat = PSI - mean(PSI),
         dPSI_logit = logit(PSI) - logit(mean(PSI)))







####  Using permutation testing ----
cat("Main tests\n")

regression_permutations <- expand_grid(event_id = unique(quantifs_filtered$event_id),
                                       method = c("lasso"),
                                       column = c("dPSI_nat"),
                                       shuffle = c(FALSE, rep(TRUE, 1000))) |>
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


cat("Main loop done. Postprocessing\n")

perm_effect_size <- regression_permutations |>
  select(event_id, shuffle, coefs_sf) |>
  unnest(coefs_sf) |>
  group_by(shuffle, event_id, transcript_id) |>
  summarize(mean_s = mean(s1),
            .groups = "drop") |>
  pivot_wider(c(event_id, transcript_id),
              names_from = "shuffle",
              values_from = "mean_s") |>
  rename(coef = `FALSE`,
         mean_null = `TRUE`) |>
  mutate(coef_effect_size = coef - mean_null)




perm_p_val <- regression_permutations |>
  select(event_id, shuffle, coefs_sf) |>
  unnest(coefs_sf) |>
  pivot_wider(c(event_id, transcript_id),
              names_from = "shuffle",
              values_from = "s1",
              values_fn = list) |>
  rename(coef = `FALSE`,
         coefs_null = `TRUE`) |>
  rowwise() |>
  mutate(p_val = sum(abs(coefs_null) >= abs(coef))/length(coefs_null)) |>
  ungroup() |>
  mutate(p_adj = p.adjust(p_val, "BH")) |>
  select(event_id, transcript_id, p_val, p_adj)



perm_res <- left_join(perm_effect_size,
                      perm_p_val,
                      by = c("event_id", "transcript_id"))

qs::qsave(perm_res, out_file)

cat("Done at ", date(),"\nsaved in ", out_file)

