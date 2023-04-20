# Compare regression on real data and on simulation


## Inits ----

library(tidyverse)
library(glmnet)

source("R/regression_functions.R")



# Predict DeltaPSI instead of PSI
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9 for rationale
logit <- function(x){
  stopifnot(all((x >= 0 & x <= 1) | is.nan(x)))
  x[x == 0] <- .01
  x[x == 1] <- .99
  log(x/(1-x))
}



## Simulated data ----


quantifs_filtered_sim <- qs::qread("data/intermediates/230413_simulation_v7/quantifs_filtered.qs")
sf_sim <- qs::qread("data/intermediates/230413_simulation_v7/sim_sf.qs")
true_coefs_sim <- qs::qread("data/intermediates/230413_simulation_v7/true_coefs.qs")



quantifs_filtered_sim <- quantifs_filtered_sim |>
  group_by(event_id) |>
  mutate(dPSI_nat = PSI - mean(PSI),
         dPSI_logit = logit(PSI) - logit(mean(PSI)))


# Make SF expression as a matrix for use in regression
mat_sf_sim <- sf_sim |>
  select(transcript_id, sample_id, TPM) |>
  mutate(TPM = log(TPM + 1)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "TPM") |>
  column_to_rownames("sample_id") |>
  as.matrix() |>
  scale()
mat_sf_sim <- mat_sf_sim[,! apply(mat_sf_sim, 2, \(col) any(is.na(col)))]


quantifs_filtered_sim
dim(mat_sf_sim)
mat_sf_sim[1:3,1:3]





## Real data ----

# read data
quantifs_filtered_real <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression_real <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs")



#preprocess data
quantifs_filtered_real <- quantifs_filtered_real |>
  group_by(event_id) |>
  mutate(dPSI_nat = PSI - mean(PSI),
         dPSI_logit = logit(PSI) - logit(mean(PSI)))

mat_sf_real <- sf_expression_real |>
  select(transcript_id, sample_id, TPM) |>
  mutate(TPM = log(TPM + 1)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "TPM") |>
  column_to_rownames("sample_id") |>
  as.matrix() |>
  scale()
mat_sf_real <- mat_sf_real[,! apply(mat_sf_real, 2, \(col) any(is.na(col)))]


quantifs_filtered_real
dim(mat_sf_real)
mat_sf_real[1:3,1:3]



## Regression ----




# First pass ----
# make a first round of regressions

first_pass_real <- expand_grid(event_id = unique(quantifs_filtered_real$event_id),
                               method = c("lasso", "adaptive_lasso"),
                               column = c("PSI", "dPSI_nat", "dPSI_logit"),
                               intercept = c(TRUE, FALSE)) |>
  # slice_sample(n = 2) |>
  mutate(res = pmap(list(event_id, method, column, intercept),
                    \(event_id, method, column, intercept) do_regression(event_id, method, column,
                                                                         shuffle = FALSE, mat_sf_real,
                                                                         quantifs_filtered_real,
                                                                         intercept),
                    .progress = TRUE),
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)




first_pass_sim <- expand_grid(event_id = unique(quantifs_filtered_sim$event_id),
                              method = c("lasso", "adaptive_lasso"),
                              column = c("PSI", "dPSI_nat", "dPSI_logit"),
                              intercept = c(TRUE, FALSE)) |>
  # slice_sample(n = 2) |>
  mutate(res = pmap(list(event_id, method, column, intercept),
                    \(event_id, method, column, intercept) do_regression(event_id, method, column,
                                                                         shuffle = FALSE, mat_sf_sim,
                                                                         quantifs_filtered_sim,
                                                                         intercept),
                    .progress = TRUE),
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)

# qs::qsave(first_pass_real, "data/intermediates/230413_simulation_v7/first_pass_real.qs")
# qs::qsave(first_pass_sim, "data/intermediates/230413_simulation_v7/first_pass_sim")


first_pass_real <- qs::qread("data/intermediates/230413_simulation_v7/first_pass_real.qs")
first_pass_sim <- qs::qread("data/intermediates/230413_simulation_v7/first_pass_sim")

table(first_pass_real$rsquare |> is.na())
table(first_pass_sim$rsquare |> is.na())

hist(first_pass_real$nb_coefs, breaks = 150)
hist(first_pass_sim$nb_coefs, breaks = 150)



hist(first_pass_real$rsquare[which(first_pass_real$nb_coefs > 1)], breaks = 50)
hist(first_pass_sim$rsquare[which(first_pass_sim$nb_coefs > 1)], breaks = 50)


bind_rows(first_pass_real |> add_column(data = "real"),
          first_pass_sim |> add_column(data = "simulated")) |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = data, y = rsquare, fill = interaction(method, column, intercept)))


bind_rows(first_pass_real |> add_column(data = "real"),
          first_pass_sim |> add_column(data = "simulated")) |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = data, y = nb_coefs, color = method, fill = interaction(column, intercept)))



# faceted version
bind_rows(first_pass_real |> add_column(data = "real") |> filter(nb_coefs > 1),
          first_pass_sim |> add_column(data = "simulated") |> filter(nb_coefs > 1)) |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = data, y = rsquare, fill = intercept)) +
  facet_grid(rows = vars(method),
             cols = vars(column))


bind_rows(first_pass_real |> add_column(data = "real"),
          first_pass_sim |> add_column(data = "simulated")) |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = data, y = nb_coefs, fill = intercept)) +
  facet_grid(rows = vars(method),
             cols = vars(column))




# Permutation tests ----

# use LASSO on deltaPSI, no intercept

# real data
permutations_real <- expand_grid(event_id = unique(quantifs_filtered_real$event_id),
                                 method = c("lasso"),
                                 column = c("dPSI_nat"),
                                 shuffle = c(rep(FALSE, 1), rep(TRUE, 100))) |>
  mutate(res = pmap(list(event_id, method, column, shuffle),
                    \(event_id, method, column, shuffle) do_regression(event_id,
                                                                       method, column,
                                                                       shuffle, mat_sf_real, quantifs_filtered_real,
                                                                       intercept = FALSE),
                    .progress = TRUE),
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)

qs::qsave(permutations_real, "data/intermediates/230413_simulation_v7/permutations_real_singleUnshuffled_100.qs")






# simulated data
permutations_sim <- expand_grid(event_id = unique(quantifs_filtered_sim$event_id),
                                 method = c("lasso"),
                                 column = c("dPSI_nat"),
                                 shuffle = c(rep(FALSE, 1), rep(TRUE, 100))) |>
  mutate(res = pmap(list(event_id, method, column, shuffle),
                    \(event_id, method, column, shuffle) do_regression(event_id,
                                                                       method, column,
                                                                       shuffle, mat_sf_sim, quantifs_filtered_sim,
                                                                       intercept = FALSE),
                    .progress = TRUE),
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)

qs::qsave(permutations_sim, "data/intermediates/230413_simulation_v7/permutations_sim_singleUnshuffled_100.qs")


#~ reformulate ----

perm_effect_size_real <- permutations_real |>
  select(event_id, shuffle, coefs_sf) |>
  unnest(coefs_sf) |>
  group_by(shuffle, event_id, transcript_id) |>
  summarize(mean_s = mean(s1),
            .groups = "drop") |>
  pivot_wider(id_cols = c(event_id, transcript_id),
              names_from = "shuffle",
              values_from = "mean_s") |>
  rename(coef = `FALSE`,
         mean_null = `TRUE`) |>
  mutate(coef_effect_size = coef - mean_null)




perm_p_val_real <- permutations_real |>
  select(event_id, shuffle, coefs_sf) |>
  unnest(coefs_sf) |>
  pivot_wider(id_cols = c(event_id, transcript_id),
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



res_real <- left_join(perm_effect_size_real,
                      perm_p_val_real,
                      by = c("event_id", "transcript_id"))





perm_effect_size_sim <- permutations_sim |>
  select(event_id, shuffle, coefs_sf) |>
  unnest(coefs_sf) |>
  group_by(shuffle, event_id, transcript_id) |>
  summarize(mean_s = mean(s1),
            .groups = "drop") |>
  pivot_wider(id_cols = c(event_id, transcript_id),
              names_from = "shuffle",
              values_from = "mean_s") |>
  rename(coef = `FALSE`,
         mean_null = `TRUE`) |>
  mutate(coef_effect_size = coef - mean_null)




perm_p_val_sim <- permutations_sim |>
  select(event_id, shuffle, coefs_sf) |>
  unnest(coefs_sf) |>
  pivot_wider(id_cols = c(event_id, transcript_id),
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



res_sim <- left_join(perm_effect_size_sim,
                      perm_p_val_sim,
                      by = c("event_id", "transcript_id"))

# Look at results

permutations_real
permutations_sim


table(is.na(permutations_real$rsquare))
hist(permutations_real$nb_coefs)

table(is.na(permutations_sim$rsquare))
hist(permutations_sim$nb_coefs)


res_real
res_sim



table(res_real$p_adj)

table(res_sim$p_adj)


res_sim |>
  filter(p_adj == 0)

xx <- res_sim |>
  left_join(true_coefs_sim,
            by = c("event_id", "transcript_id"))


xx |>
  anti_join(true_coefs_sim, by = c("event_id", "transcript_id")) |> pull(p_adj) |> table()
true_coefs_sim |>
  anti_join(xx, by = c("event_id", "transcript_id")) |>
  pull(true_coef) -> yy

table(abs(yy) > 1)

xx |>
  filter(!is.na(true_coef)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(true_coef, y = coef_effect_size, shape = contribution, color = p_adj == 0))
  


