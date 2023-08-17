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



## Simulated data v11 ----


quantifs_filtered_sim <- qs::qread("data/intermediates/230517_simulation_v12/quantifs_filtered.qs")
sf_sim <- qs::qread("data/intermediates/230517_simulation_v12/sim_sf.qs")
true_coefs_sim <- qs::qread("data/intermediates/230517_simulation_v12/true_coefs.qs")



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



## Simulated data v1 ----


quantifs_filtered_sim1 <- qs::qread("data/intermediates/simultation/sim_quantifs.qs")
sf_sim1 <- qs::qread("data/intermediates/simultation/sim_sf.qs")
true_coefs_sim1 <- qs::qread("data/intermediates/simultation/true_coefs.qs")



quantifs_filtered_sim1 <- quantifs_filtered_sim1 |>
  group_by(event_id) |>
  mutate(dPSI_nat = PSI - mean(PSI),
         dPSI_logit = logit(PSI) - logit(mean(PSI)))


# Make SF expression as a matrix for use in regression
mat_sf_sim1 <- sf_sim1 |>
  select(transcript_id, sample_id, TPM) |>
  mutate(TPM = log(TPM + 1)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "TPM") |>
  column_to_rownames("sample_id") |>
  as.matrix() |>
  scale()
mat_sf_sim1 <- mat_sf_sim1[,! apply(mat_sf_sim1, 2, \(col) any(is.na(col)))]


quantifs_filtered_sim1
dim(mat_sf_sim1)
mat_sf_sim1[1:3,1:3]





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
                               method = c("adaptive_lasso"),
                               column = c("dPSI_nat"),
                               intercept = c(FALSE)) |>
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
                              method = c("adaptive_lasso"),
                              column = c("dPSI_nat"),
                              intercept = c(FALSE)) |>
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

qs::qsave(first_pass_real, "data/intermediates/230517_simulation_v11/first_pass_real.qs")
qs::qsave(first_pass_sim, "data/intermediates/230517_simulation_v11/first_pass_sim")

#~ Check results ----
# first_pass_real <- qs::qread("data/intermediates/230517_simulation_v11/first_pass_real.qs")
# first_pass_sim <- qs::qread("data/intermediates/230517_simulation_v11/first_pass_sim")

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
             cols = vars(column)) +
  geom_hline(aes(yintercept = 40), color = 'grey', linetype = 'dashed')






# Check fit on training set ----

# Use adaptive lasso on deltaPSI no intercept

#~ Real data ----

# (my_ev <- sample(unique(quantifs_filtered_real$event_id), 1))
rsqu_real <- double(length(unique(quantifs_filtered_real$event_id))) |>
  setNames(unique(quantifs_filtered_real$event_id))
for(my_ev in unique(quantifs_filtered_real$event_id)){
  # get y data
  y <- quantifs_filtered_real[quantifs_filtered_real$event_id == my_ev, c("sample_id", "dPSI_nat")] |>
    column_to_rownames("sample_id") |>
    filter(!is.na(.data[["dPSI_nat"]])) |>
    as.matrix()
  
  # get x data
  x <- mat_sf_real[rownames(y),]
  
  n <- nrow(x)
  train <- sample(n, round(.7*n))
  
  
  
  fit <- adaptive_lasso(x[train,], y[train], nfolds = 20, intercept = FALSE)
  
  # Look in train data
  y_predicted <- predict(fit, newx = x[train,], s = "lambda.min")
  
  # plot(y[train,], y_predicted, xlab = "PSI (measured)", ylab = "PSI (from fit)")
  
  rsqu_real[[my_ev]] <- summary(lm(y_predicted ~ y[train,]))$adj.r.squared
}


#~ Simulated data v12 ----
# (my_ev <- sample(unique(quantifs_filtered_sim$event_id), 1))
rsqu_sim <- double(length(unique(quantifs_filtered_sim$event_id))) |>
  setNames(unique(quantifs_filtered_sim$event_id))
for(my_ev in unique(quantifs_filtered_sim$event_id)){
  # get y data
  y <- quantifs_filtered_sim[quantifs_filtered_sim$event_id == my_ev, c("sample_id", "dPSI_nat")] |>
    column_to_rownames("sample_id") |>
    filter(!is.na(.data[["dPSI_nat"]])) |>
    as.matrix()
    # get x data
  x <- mat_sf_sim[rownames(y),]
    n <- nrow(x)
  train <- sample(n, round(.7*n))
    fit <- adaptive_lasso(x[train,], y[train], nfolds = 20, intercept = FALSE)
  
  # Look in train data
  y_predicted <- predict(fit, newx = x[train,], s = "lambda.min")
  
  # plot(y[train,], y_predicted, xlab = "PSI (simulated)", ylab = "PSI (from fit)")

  rsqu_sim[[my_ev]] <- summary(lm(y_predicted ~ y[train,]))$adj.r.squared
}

table(is.nan(rsqu_real))
table(is.nan(rsqu_sim))

boxplot(c(rsqu_real, rsqu_sim) ~ rep(c("real", "simulated"),
                                     times = c(length(rsqu_real), length(rsqu_sim))),
        xlab = NULL, ylab = "R² (on training data)")
abline(h = .5, col = 'grey', lty = 'dashed')




# Look at correlation ----
# Question: is correlation higher for true causal SFs?

#~ sim ----

# check that all needed sample_id exist in mat_sf
# for(my_ev in unique(quantifs_filtered_sim$event_id)){
#   y <- quantifs_filtered_sim |>
#     filter(event_id == my_ev) |>
#     pull(sample_id)
#   stopifnot(all(y %in% rownames(mat_sf_sim)))
# }

# for each event, for each SF, across samples
correlations_sim <- map(unique(quantifs_filtered_sim$event_id),
    \(.ev) {
      y <- quantifs_filtered_sim |>
        filter(event_id == .ev)
      
      cor(mat_sf_sim[y$sample_id, ],
          y$PSI) |>
        as.data.frame() |>
        rownames_to_column("transcript_id") |>
        add_column(event_id = .ev)
      
      }) |>
  list_rbind() |>
  as_tibble() |>
  rename(R = V1) |>
  left_join(true_coefs_sim,
            by = c("event_id", "transcript_id")) |>
  mutate(relevant = true_coef != 0)


correlations_sim |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = abs(R), fill = relevant, color = relevant),
               alpha = .5, size = 1) +
  # scale_fill_manual(values = c('#5182AF', '#C04045')) +
  # scale_color_manual(values = c('#5182AF', '#C04045'))
  theme(legend.position = c(0.8,.5))

correlations_sim |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = abs(R), y = runif(nrow(correlations_sim), 0, 100), color = relevant)) +
  scale_color_manual(values = c('#5182AF', '#C04045'))

correlations_real <- map(unique(quantifs_filtered_real$event_id),
                        \(.ev) {
                          y <- quantifs_filtered_real |>
                            filter(event_id == .ev)
                          
                          cor(mat_sf_real[y$sample_id, ],
                              y$PSI) |>
                            as.data.frame() |>
                            rownames_to_column("transcript_id") |>
                            add_column(event_id = .ev)
                          
                        }) |>
  list_rbind() |>
  as_tibble() |>
  rename(R = V1)

correlations_sim |>
  add_column(data = "simulated") |>
  bind_rows(correlations_real |> add_column(data = "real")) |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = abs(R), fill = interaction(data, relevant)),
               alpha = .5) +
  theme(legend.position = c(0.8,.5))



# Permutation tests ----

# use adaptive LASSO on deltaPSI, no intercept

# real data
permutations_real <- expand_grid(event_id = unique(quantifs_filtered_real$event_id),
                                 method = c("adaptive_lasso"),
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

# qs::qsave(permutations_real, "data/intermediates/230512_simulation_v10/permutations_real_singleUnshuffled_100.qs")






# simulated data v7
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

# qs::qsave(permutations_sim, "data/intermediates/230512_simulation_v10/permutations_sim_singleUnshuffled_100.qs")


# simulated data v1
permutations_sim1 <- expand_grid(event_id = unique(quantifs_filtered_sim1$event_id),
                                method = c("lasso"),
                                column = c("dPSI_nat"),
                                shuffle = c(rep(FALSE, 1), rep(TRUE, 100))) |>
  mutate(res = pmap(list(event_id, method, column, shuffle),
                    \(event_id, method, column, shuffle) do_regression(event_id,
                                                                       method, column,
                                                                       shuffle, mat_sf_sim1, quantifs_filtered_sim1,
                                                                       intercept = FALSE),
                    .progress = TRUE),
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)

# qs::qsave(permutations_sim1, "data/intermediates/230512_simulation_v10/permutations_sim1_singleUnshuffled_100.qs")




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



perm_effect_size_sim1 <- permutations_sim1 |>
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




perm_p_val_sim1 <- permutations_sim1 |>
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



res_sim1 <- left_join(perm_effect_size_sim1,
                     perm_p_val_sim1,
                     by = c("event_id", "transcript_id"))



# Look at results

permutations_real
permutations_sim
permutations_sim1

table(is.na(permutations_real$rsquare))
hist(permutations_real$nb_coefs)

table(is.na(permutations_sim$rsquare))
hist(permutations_sim$nb_coefs)


table(is.na(permutations_sim1$rsquare))
hist(permutations_sim1$nb_coefs)

res_real
res_sim
res_sim1


table(res_real$p_adj)

table(res_sim$p_adj)
table(res_sim1$p_adj)


res_sim |>
  filter(p_adj == 0)

xx <- res_sim1 |>
  left_join(true_coefs_sim1,
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
  geom_point(aes(true_coef, y = coef_effect_size, color = p_adj == 0))
  




# Network ----
library(DPM)

#~ sim v11 ----
y_sim <- quantifs_filtered_sim |>
  select(event_id, sample_id, dPSI_nat) |>
  pivot_wider(id_cols = event_id, names_from = sample_id, values_from = dPSI_nat) |>
  column_to_rownames("event_id") |>
  as.matrix()
dim(y_sim)
y_sim[1:3,1:3]

# get x data
dat_sim <- cbind(mat_sf_sim,t(y_sim))


dim(dat_sim)
dat_sim[1:3,1:3]

table(is.na(dat_sim))
dat_sim2 <- impute::impute.knn(dat_sim)$data
dpm_sim <- DPM::reg.dpm(dat_sim2)

# qs::qsave(dpm_sim, "data/intermediates/230607_simulation_v16/dpm_sim.qs")
dim(dpm_sim)
dpm_sim[1:4,1:4]
diag(dpm_sim) <- 0
image(dpm_sim)
rownames(dpm_sim) <- colnames(dat_sim)
colnames(dpm_sim) <- colnames(dat_sim)


nb_sf <- ncol(mat_sf_sim)
nb_psi <- nrow(y_sim)
stopifnot(nrow(dpm_sim) == nb_sf + nb_psi)

annot <- data.frame(row.names = rownames(dpm_sim),
                    type = rep(c("SF", "PSI"),
                               times = c(nb_sf, nb_psi)))

crosscorr_sim <- dpm_sim[1:nb_sf, (nb_sf + 1):(nb_sf + nb_psi)]

dim(crosscorr_sim)

all.equal(t(crosscorr_sim),
          dpm_sim[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf])

image(crosscorr_sim)



pheatmap::pheatmap(crosscorr_sim,
                   show_rownames = FALSE,
                   show_colnames = FALSE)



#~ real ----
y_real <- quantifs_filtered_real |>
  select(event_id, sample_id, dPSI_nat) |>
  pivot_wider(id_cols = event_id, names_from = sample_id, values_from = dPSI_nat) |>
  column_to_rownames("event_id") |>
  as.matrix()
dim(y_real)
y_real[1:3,1:3]

# get x data
dat_real <- cbind(mat_sf_real,t(y_real))


dim(dat_real)
dat_real[1:3,1:3]

table(is.na(dat_real))
dat_real2 <- impute::impute.knn(dat_real)$data
dpm_real <- DPM::reg.dpm(dat_real2)

# qs::qsave(dpm_real, "data/intermediates/230607_simulation_v16/dpm_real.qs")
dpm_real <- qs::qread("data/intermediates/230607_simulation_v16/dpm_real.qs")

dim(dpm_real)
dpm_real[1:4,1:4]
diag(dpm_real) <- 0
image(dpm_real)
rownames(dpm_real) <- colnames(dat_real)
colnames(dpm_real) <- colnames(dat_real)


nb_sf <- ncol(mat_sf_real)
nb_psi <- nrow(y_real)
stopifnot(nrow(dpm_real) == nb_sf + nb_psi)

crosscorr_real <- dpm_real[1:nb_sf, (nb_sf + 1):(nb_sf + nb_psi)]

dim(crosscorr_real)

all.equal(t(crosscorr_real),
          dpm_real[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf])

image(crosscorr_real)


pheatmap::pheatmap(crosscorr_real,
                   show_rownames = FALSE,
                   show_colnames = FALSE)




dpm_links <- DPM::get_link(dpm_real)

dpm_links |>
  mutate(type1 = if_else(Node1 %in% colnames(mat_sf_real),
                         "SF", "PSI"),
         type2 = if_else(Node2 %in% colnames(mat_sf_real),
                         "SF", "PSI"),
         bipartite = type1 != type2) |>
  pull(bipartite) |> table()

xx <- dpm_links |>
  as_tibble() |>
  mutate(type1 = if_else(Node1 %in% colnames(mat_sf_real),
                         "SF", "PSI"),
         type2 = if_else(Node2 %in% colnames(mat_sf_real),
                         "SF", "PSI"),
         bipartite = type1 != type2) |>
  filter(bipartite) |>
  arrange(desc(score))

hist(log10(xx$score), breaks = 50); abline(v = log10(2e-2), col = 'darkred')

xx2 <- xx[xx$score > 2e-2,]
qgraph::qgraph(xx2)

regdpmdat <-tibble(score = abs( dpm_real[lower.tri(dpm_real)]), name='DPM_real')

ggplot(regdpmdat, aes(x=score,color=name))+ theme_bw()  + geom_density(size=1.5)+
  scale_color_manual(values=c("DPM"="SteelBlue3","DPM_real"="dodgerblue4"))

regdpmtr <- kmeans_links(dpm_real)$threshold


library(igraph)
library(Matrix)

sp_crosscorr <- Matrix(crosscorr_real, sparse = TRUE)
sp_crosscorr[sp_crosscorr < 2e-2] <- 0

names_sf <- colnames(mat_sf_real)
names_psi <- rownames(y_real)

adjm <- rbind(
  cbind(Matrix(0, nrow = nb_sf, ncol = nb_sf, dimnames = list(names_sf, names_sf)),
        sp_crosscorr),
  cbind(Matrix(0, nrow = nb_psi, ncol = nb_sf, dimnames = list(names_psi, names_sf)),
        Matrix(0, nrow = nb_psi, ncol = nb_psi, dimnames = list(names_psi, names_psi)))
)

image(adjm)
image(adjm[1:nb_sf, (nb_sf+1):(nb_sf+nb_psi)])


gr <- graph_from_adjacency_matrix(adjm,
                                  mode = "directed",
                                  weighted = TRUE) |>
  set_vertex_attr("type",
                  value = rep(c("sf", "psi"), times = c(nb_sf, nb_psi)))
gr

length(V(gr))
nb_sf+nb_psi

length(E(gr))
sum(adjm > 0)

head(E(gr))
tail(E(gr))

RCy3::cytoscapePing()
RCy3::createNetworkFromIgraph(gr)




plot(gr,
     layout = layout_as_bipartite,
     vertex.color = V(gr)$type,
     label.cex = .001)










