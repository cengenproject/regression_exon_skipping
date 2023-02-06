# Simulation to check whether the permutation test gives results that make sense



# ----Simulations ----

library(tidyverse)

# read filtered objects to generate simulations from
quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs")

logistic <- function(x, L = 1, k = 1, x0 = 0){
  L/(1+exp(-k*(x-x0)))
}


# simulate single
sim_quantifs <- quantifs_filtered |>
  select(event_id, sample_id)
sim_sf <- sf_expression |>
  select(transcript_id, sample_id) |>
  mutate(TPM = sample(sf_expression$TPM))

n_tx <- length(unique(sim_sf$transcript_id))
true_coefs <- expand_grid(event_id = unique(sim_quantifs$event_id),
                          transcript_id = unique(sim_sf$transcript_id)) |>
  group_by(event_id) |>
  nest() |>
  mutate(data = map(data, ~ add_column(.x, true_coef = sample(c(rep(0, n_tx - 5), rnorm(5)))))) |>
  unnest(data) |> ungroup()


sim_quantifs <- sim_quantifs |>
  left_join(sim_sf, by = "sample_id") |>
  left_join(true_coefs, by = c("event_id", "transcript_id")) |>
  group_by(event_id, sample_id) |>
  summarize(PSI = sum(true_coef*TPM),
            .groups = 'drop') |>
  mutate(PSI = PSI + rnorm(length(PSI), mean = 0, sd = .2),
         PSI = logistic(PSI, k = .05))

qs::qsave(sim_quantifs, "data/intermediates/simultation/sim_quantifs.qs")
qs::qsave(sim_sf, "data/intermediates/simultation/sim_sf.qs")
qs::qsave(true_coefs, "data/intermediates/simultation/true_coefs.qs")




# Check simulated PSI
(my_ev <- sample(events_to_keep, 1))
bind_rows(
  quantifs_filtered |>
    filter(event_id == my_ev) |>
    select(ends_with("_id"), PSI) |>
    add_column(type = "measured"),
  sim_quantifs |>
    filter(event_id == my_ev) |>
    select(ends_with("_id"), PSI) |>
    mutate(neuron_id = str_match(sample_id, "^([A-Z0-9]+)r[0-9]+$")[,2]) |>
    add_column(type = "simulated")) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~type) +
  geom_col(aes(x = sample_id, y = PSI, fill = neuron_id, group = neuron_id),
           position = "dodge", show.legend = FALSE, color = 'white') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
  hues::scale_fill_iwanthue()

hist(quantifs_filtered$PSI, breaks = 150)
hist(sim_quantifs$PSI, breaks = 150)

bind_rows(
  quantifs_filtered |>
    select(ends_with("_id"), PSI) |>
    add_column(type = "measured"),
  sim_quantifs |>
    select(ends_with("_id"), PSI) |>
    mutate(neuron_id = str_match(sample_id, "^([A-Z0-9]+)r[0-9]+$")[,2]) |>
    add_column(type = "simulated")) |>
  ggplot() +
  theme_classic() +
  geom_freqpoly(aes(x = PSI, color = type), bins = 100)



# Simulate 100 datasets


sim_replicated <- replicate(100,
                            {sim_quantifs <- quantifs_filtered |>
                              select(event_id, sample_id)
                            
                            sim_sf <- sf_expression |>
                              select(transcript_id, sample_id) |>
                              mutate(TPM = sample(sf_expression$TPM))
                            
                            n_tx <- length(unique(sim_sf$transcript_id))
                            
                            true_coefs <- expand_grid(event_id = unique(sim_quantifs$event_id),
                                                      transcript_id = unique(sim_sf$transcript_id)) |>
                              group_by(event_id) |>
                              nest() |>
                              mutate(data = map(data, ~ add_column(.x, true_coef = sample(c(rep(0, n_tx - 5), rnorm(5)))))) |>
                              unnest(data) |> ungroup()
                            
                            
                            sim_quantifs <- sim_quantifs |>
                              left_join(sim_sf, by = "sample_id") |>
                              left_join(true_coefs, by = c("event_id", "transcript_id")) |>
                              group_by(event_id, sample_id) |>
                              summarize(PSI = sum(true_coef*TPM),
                                        .groups = 'drop') |>
                              mutate(PSI = PSI + rnorm(length(PSI), mean = 0, sd = .2),
                                     PSI = logistic(PSI, k = .05))
                            list(sim_quantifs = sim_quantifs, sim_sf = sim_sf, true_coefs = true_coefs)
                            },
                            simplify = FALSE)


# qs::qsave(sim_replicated, "data/intermediates/simultation/230206_rep_simulations.qs")

bind_rows(
  quantifs_filtered |>
    select(ends_with("_id"), PSI) |>
    add_column(type = "measured"),
  sim_replicated[[2]][[1]] |>
    select(ends_with("_id"), PSI) |>
    mutate(neuron_id = str_match(sample_id, "^([A-Z0-9]+)r[0-9]+$")[,2]) |>
    add_column(type = "simulated")) |>
  ggplot() +
  theme_classic() +
  geom_freqpoly(aes(x = PSI, color = type), bins = 100)

# Export 2 replicates for permutation tests
# qs::qsave(sim_replicated[[45]]$sim_quantifs, "data/intermediates/simultation/230206_quantifs_simulation1.qs")
# qs::qsave(sim_replicated[[45]]$sim_sf, "data/intermediates/simultation/230206_sf_simulation1.qs")
# 
# qs::qsave(sim_replicated[[10]]$sim_quantifs, "data/intermediates/simultation/230206_quantifs_simulation2.qs")
# qs::qsave(sim_replicated[[10]]$sim_sf, "data/intermediates/simultation/230206_sf_simulation2.qs")





# ---- Regression ----

library(tidyverse)
library(glmnet)

source("R/regression_functions.R")


# Read data ----
sim_replicated <- qs::qread("data/intermediates/simultation/230206_rep_simulations.qs")

sim_quantifs <- list_transpose(sim_replicated)[["sim_quantifs"]]
sim_sf <- list_transpose(sim_replicated)[["sim_sf"]]
sim_true_coefs <- list_transpose(sim_replicated)[["true_coefs"]]



# Predict DeltaPSI instead of PSI
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9 for rationale

logit <- function(x){
  stopifnot(all(x >= 0 & x <= 1))
  x[x == 0] <- .01
  x[x == 1] <- .99
  log(x/(1-x))
}

sim_quantifs <- map(sim_quantifs,
                    ~ .x |>
                      group_by(event_id) |>
                      mutate(dPSI_nat = PSI - mean(PSI),
                             dPSI_logit = logit(PSI) - logit(mean(PSI))))


events_to_keep <- unique(sim_quantifs[[1]]$event_id)





# Make SF expression as a matrix for use in regression

sim_mat_sf <- map(sim_sf,
              ~ {mat_sf_expression <-  .x|>
                  select(transcript_id, sample_id, TPM) |>
                  mutate(TPM = log(TPM + 1)) |>
                  pivot_wider(id_cols = sample_id,
                              names_from = "transcript_id",
                              values_from = "TPM") |>
                  column_to_rownames("sample_id") |>
                  as.matrix() |>
                  scale()
                mat_sf_expression[,! apply(mat_sf_expression, 2, \(col) any(is.na(col)))]})







###

# Play with examples ----
# SE_136, SE_303 are good. 1060 particularly
(my_ev <- sample(events_to_keep, 1))
# my_ev <- "SE_1046"


sim_quantifs[[37]] |>
  filter(event_id == my_ev) |>
  pivot_longer(-ends_with("_id"),
               names_to = "type",
               values_to = "PSI") |>
  ggplot() +
  theme_classic() +
  facet_wrap(~type, scales = "free_y") +
  geom_col(aes(x = sample_id, y = PSI, fill = sample_id, group = sample_id),
           position = "dodge", show.legend = FALSE, color = 'white') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
  hues::scale_fill_iwanthue()















# For single replicate ----
quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_rep_sim_quantifs.qs")
sf_expression <- qs::qread("data/intermediates/simultation/sim_sf.qs")




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








###
events_to_keep <- unique(quantifs_filtered$event_id)

# Play with examples ----
# SE_136, SE_303 are good. 1060 particularly
(my_ev <- sample(events_to_keep, 1))
# my_ev <- "SE_1046"


quantifs_filtered |>
  filter(event_id == my_ev) |>
  pivot_longer(-ends_with("_id"),
               names_to = "type",
               values_to = "PSI") |>
  ggplot() +
  theme_classic() +
  facet_wrap(~type, scales = "free_y") +
  geom_col(aes(x = sample_id, y = PSI, fill = sample_id, group = sample_id),
           position = "dodge", show.legend = FALSE, color = 'white') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
  hues::scale_fill_iwanthue()






# Sparse LASSO


reg_lasso <- regression_wrapper(my_ev = my_ev, regression_method = "lasso", column = "PSI",
                                shuffle = FALSE, mat_sf_expression = mat_sf_expression, quants = quantifs_filtered)
reg_adalasso <- regression_wrapper(my_ev = my_ev, regression_method = "adaptive_lasso", column = "PSI",
                                   shuffle = FALSE, mat_sf_expression = mat_sf_expression, quants = quantifs_filtered)
reg_lasso_nat <- regression_wrapper(my_ev = my_ev, regression_method = "lasso", column = "dPSI_nat",
                                    shuffle = FALSE, mat_sf_expression = mat_sf_expression, quants = quantifs_filtered)
reg_adalasso_nat <- regression_wrapper(my_ev = my_ev, regression_method = "adaptive_lasso", column = "dPSI_nat",
                                       shuffle = FALSE, mat_sf_expression = mat_sf_expression, quants = quantifs_filtered)
reg_lasso_logit <- regression_wrapper(my_ev = my_ev, regression_method = "lasso", column = "dPSI_logit",
                                      shuffle = FALSE, mat_sf_expression = mat_sf_expression, quants = quantifs_filtered)
reg_adalasso_logit <- regression_wrapper(my_ev = my_ev, regression_method = "adaptive_lasso", column = "dPSI_logit",
                                         shuffle = FALSE, mat_sf_expression = mat_sf_expression, quants = quantifs_filtered)




# On test set
reg_lasso$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")

reg_adalasso$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")

reg_lasso_nat$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")

reg_adalasso_nat$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")

reg_lasso_logit$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")

reg_adalasso_logit$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")



# coefficients
reg_lasso$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = transcript_id, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


reg_adalasso$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = transcript_id, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 

reg_lasso_nat$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = transcript_id, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


reg_adalasso_nat$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = transcript_id, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 

reg_lasso_logit$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = transcript_id, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


reg_adalasso_logit$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = transcript_id, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


list(lasso    = reg_lasso$coefs_sf |>  filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     adalasso = reg_adalasso$coefs_sf |> filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     lasso_nat  = reg_lasso_nat$coefs_sf |>  filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     adalasso_nat = reg_adalasso_nat$coefs_sf |> filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     lasso_logit  = reg_lasso_logit$coefs_sf |>  filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     adalasso_logit = reg_adalasso_logit$coefs_sf |> filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)")) |>
  eulerr::euler() |>
  plot(quantities = TRUE)



# Compare to true
true_coefs <- qs::qread("data/intermediates/simultation/true_coefs.qs")

comp_to_true <- list(lasso    = reg_lasso$coefs_sf |> 
                       add_column(method = "lasso"),
                     adalasso = reg_adalasso$coefs_sf |> 
                       add_column(method = "adalasso"),
                     lasso_nat  = reg_lasso_nat$coefs_sf |> 
                       add_column(method = "lasso_nat"),
                     adalasso_nat = reg_adalasso_nat$coefs_sf |> 
                       add_column(method = "adalasso_nat"),
                     lasso_logit  = reg_lasso_logit$coefs_sf |>
                       add_column(method = "lasso_logit"),
                     adalasso_logit = reg_adalasso_logit$coefs_sf |>
                       add_column(method = "adalasso_logit")) |>
  list_rbind() |>
  add_column(event_id = my_ev) |>
  left_join(true_coefs, by = c("event_id", "transcript_id"))


comp_to_true |>
  filter(method == "lasso_logit") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = true_coef, y = s1)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  coord_equal()





