
## ARCHIVE


# Simulation to check whether the permutation test gives results that make sense

# While the simulated data appears good enough at first, it has deviations from
# the real data (e.g. the PSI distribution has a peak at 0, and a small peak at 0.5),
# so replacing with new script, see `simulation2.R`


# Model in this script:
# Keep measured TPMs, just randomize them
# For each event, select 5 transcripts non-zero, all other coefficients 0
# For these 5 “true non-0” transcripts, coefficient from Norm(mu=0, sd=1)
# Compute simulated
# PSI = logistic( sum(coef * TPM) + eps )
# where eps from Norm(mu = 0, sd = 0.2), logistic with k=0.05
# Create 100 simulated datasets
# 
# Note: in model, use scaled-centered log-TPM





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




#~ Look at single replicate ----

#~ lambda selection ----
# Check lambda for single replicate, single event
(my_ev <- sample(events_to_keep, 1))
(my_rep <- sample(100, 1))
quantifs_filtered <- sim_quantifs[[my_rep]]
mat_sf_expression <- sim_mat_sf[[my_rep]]
y <- quantifs_filtered[quantifs_filtered$event_id == my_ev, c("sample_id", "PSI")] |>
  column_to_rownames("sample_id") |>
  filter(!is.na("PSI")) |>
  as.matrix()

# get x data
x <- mat_sf_expression[rownames(y),]

n <- nrow(x)
train <- sample(n, round(.7*n))

fit <- cv.glmnet(x = x[train,], y = y[train,], type.measure = "mse", nfolds = 20, alpha = 1)
plot(fit)

log(fit$lambda.min)
log(fit$lambda.1se)



# Sparse LASSO

# sim1
quantifs_filtered <- sim_quantifs[[1]]
mat_sf_expression <- sim_mat_sf[[1]]

reg_lasso1 <- regression_wrapper(my_ev = my_ev, regression_method = "lasso", column = "dPSI_nat",
                                 shuffle = FALSE, mat_sf_expression = mat_sf_expression, quants = quantifs_filtered)

reg_lasso1$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  coord_fixed()

reg_lasso1$coefs_sf |>
  filter(s1 != 0) |>
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
true_coefs1 <- sim_true_coefs[[1]] |> filter(event_id == my_ev)


# sim 2
quantifs_filtered <- sim_quantifs[[2]]
mat_sf_expression <- sim_mat_sf[[2]]
reg_lasso2 <- regression_wrapper(my_ev = my_ev, regression_method = "lasso", column = "dPSI_nat",
                                 shuffle = FALSE, mat_sf_expression = mat_sf_expression, quants = quantifs_filtered)

reg_lasso2$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  coord_fixed()
reg_lasso2$coefs_sf |>
  filter(s1 != 0) |>
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
true_coefs2 <- sim_true_coefs[[2]] |> filter(event_id == my_ev)






comp_to_true <- list(sim1    = reg_lasso1$coefs_sf |> 
                       add_column(sim_rep = "sim1"),
                     sim2 = reg_lasso2$coefs_sf |> 
                       add_column(sim_rep = "sim2")) |>
  list_rbind() |>
  add_column(event_id = my_ev) |>
  left_join(bind_rows(true_coefs1 |> add_column(sim_rep = "sim1"),
                      true_coefs2 |> add_column(sim_rep = "sim2")),
            by = c("event_id", "transcript_id", "sim_rep"))


# rsquares
comp_to_true |>
  filter(transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = true_coef, y = s1)) +
  # geom_abline(aes(slope = 1, intercept = 0)) +
  # coord_equal() +
  facet_wrap(~sim_rep) +
  ylab("Computed coef")


map_dbl(unique(comp_to_true$sim_rep),
        ~ lm(s1 ~ true_coef, data = comp_to_true |>
               filter(transcript_id != "(Intercept)", sim_rep == .x)) |>
          summary() |>
          {\(x) x$adj.r.squared}())

# TPR FDR
comp_to_true |>
  filter(transcript_id != "(Intercept)") |>
  mutate(s1 = s1 != 0,
         true_coef = true_coef != 0) |>
  group_by(sim_rep, event_id) |>
  summarize(TP = sum(s1 & true_coef),
            FP = sum(s1 & !true_coef),
            FN = sum(!s1 & true_coef),
            TN = sum(!s1 & !true_coef),
            .groups = 'drop') |>
  mutate(TPR = TP/(TP+FN),
         FPR = FP/(FP+TN),
         FDR = FP/(FP+TP))



#~  all replicates ----

reg_lasso_coefs <- map(sample(100, 50),
                 ~ regression_wrapper(my_ev = my_ev, regression_method = "lasso", column = "PSI",
                                 shuffle = FALSE,
                                 mat_sf_expression = sim_mat_sf[[.x]],
                                 quants = sim_quantifs[[.x]]) |>
                   pluck("coefs_sf") |> 
                   add_column(sim_rep = paste0("sim_", .x)),
                 .progress = TRUE) |>
  list_rbind()


true_coefs <- imap(sim_true_coefs,
                      ~ .x |>
                        filter(event_id == my_ev) |>
                        add_column(sim_rep = paste0("sim_", .y))) |>
  list_rbind()




comp_to_true <- left_join(reg_lasso_coefs |>
                            filter(transcript_id != "(Intercept)"),
                          true_coefs |>
                            select(-event_id),
            by = c("transcript_id", "sim_rep"))


comp_to_true |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = true_coef, y = s1)) +
  # geom_abline(aes(slope = 1, intercept = 0)) +
  # coord_equal() +
  # facet_wrap(~sim_rep) +
  ylab("Computed coef")




all_rsquared <- map_dbl(unique(comp_to_true$sim_rep),
        ~ lm(s1 ~ true_coef, data = comp_to_true |>
               filter(transcript_id != "(Intercept)", sim_rep == .x)) |>
          summary() |>
          {\(x) x$adj.r.squared}())

hist(all_rsquared, breaks = 50)

# TPR FDR
comp_to_true |>
  filter(transcript_id != "(Intercept)") |>
  mutate(s1 = s1 != 0,
         true_coef = true_coef != 0) |>
  group_by(sim_rep) |>
  summarize(TP = sum(s1 & true_coef),
            FP = sum(s1 & !true_coef),
            FN = sum(!s1 & true_coef),
            TN = sum(!s1 & !true_coef),
            .groups = 'drop') |>
  mutate(TPR = TP/(TP+FN),
         FPR = FP/(FP+TN),
         FDR = FP/(FP+TP))


#~  all replicates, all events ----

get_rsquared <- function(.method, .unit){
  map(sample(events_to_keep, 10),
      \(my_ev){
        reg_lasso_coefs <- map(sample(100, 10),
                               ~ regression_wrapper(my_ev = my_ev, regression_method = .method,
                                                    column = .unit,
                                                    shuffle = FALSE,
                                                    mat_sf_expression = sim_mat_sf[[.x]],
                                                    quants = sim_quantifs[[.x]]) |>
                                 pluck("coefs_sf") |> 
                                 add_column(sim_rep = paste0("sim_", .x))) |>
          list_rbind()
        
        
        true_coefs <- imap(sim_true_coefs,
                           ~ .x |>
                             filter(event_id == my_ev) |>
                             add_column(sim_rep = paste0("sim_", .y))) |>
          list_rbind()
        
        
        
        
        comp_to_true <- left_join(reg_lasso_coefs |>
                                    filter(transcript_id != "(Intercept)"),
                                  true_coefs |>
                                    select(-event_id),
                                  by = c("transcript_id", "sim_rep"))
        
        map_dbl(unique(comp_to_true$sim_rep),
                ~ lm(s1 ~ true_coef, data = comp_to_true |>
                       filter(transcript_id != "(Intercept)", sim_rep == .x)) |>
                  summary() |>
                  {\(x) x$adj.r.squared}())
      })
}


all_rsquared <- tibble(regression_method = rep(c("lasso", "adaptive_lasso"), each = 3),
                      unit = rep(c("PSI","dPSI_nat","dPSI_logit"), times = 2)) |>
  mutate(res = map2(regression_method, unit, get_rsquared,
                    .progress = TRUE))

all_rsquared |>
  mutate(res = map(res, unlist))  |>
  unnest(res) |>
  rename(rsquared = res) |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = regression_method, fill = unit, y = rsquared)) +
  scale_y_continuous(limits = c(0,1))


# this was all using Rsquared as output, what if using the TPR and FDR


get_TPR_FDR <- function(.method, .unit){
  map(sample(events_to_keep, 10),
      \(my_ev){
        reg_lasso_coefs <- map(sample(100, 10),
                               ~ regression_wrapper(my_ev = my_ev, regression_method = .method,
                                                    column = .unit,
                                                    shuffle = FALSE,
                                                    mat_sf_expression = sim_mat_sf[[.x]],
                                                    quants = sim_quantifs[[.x]]) |>
                                 pluck("coefs_sf") |> 
                                 add_column(sim_rep = paste0("sim_", .x))) |>
          list_rbind()
        
        
        true_coefs <- imap(sim_true_coefs,
                           ~ .x |>
                             filter(event_id == my_ev) |>
                             add_column(sim_rep = paste0("sim_", .y))) |>
          list_rbind()
        
        
        
        
        comp_to_true <- left_join(reg_lasso_coefs |>
                                    filter(transcript_id != "(Intercept)"),
                                  true_coefs |>
                                    select(-event_id),
                                  by = c("transcript_id", "sim_rep"))
        
        # TPR FDR
        comp_to_true |>
          filter(transcript_id != "(Intercept)") |>
          mutate(s1 = s1 != 0,
                 true_coef = true_coef != 0) |>
          group_by(sim_rep) |>
          summarize(TP = sum(s1 & true_coef),
                    FP = sum(s1 & !true_coef),
                    FN = sum(!s1 & true_coef),
                    TN = sum(!s1 & !true_coef),
                    .groups = 'drop') |>
          mutate(TPR = TP/(TP+FN),
                 FPR = FP/(FP+TN),
                 FDR = FP/(FP+TP))
      })
}


all_TPR_FDR <- tibble(regression_method = rep(c("lasso", "adaptive_lasso"), each = 3),
                       unit = rep(c("PSI","dPSI_nat","dPSI_logit"), times = 2)) |>
  mutate(res = map2(regression_method, unit, get_TPR_FDR,
                    .progress = TRUE))


all_TPR_FDR |>
  mutate(res = map(res, list_rbind)) |>
  unnest(res) |>
  select(regression_method, unit, TPR, FPR, FDR) |>
  pivot_longer(-c(regression_method, unit),
               names_to = "metric") |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = regression_method, fill = unit, y = value)) +
  facet_wrap(~metric) +
  xlab(NULL)+ylab(NULL)




#~ Effect of intercept ----

# Check lambda for single replicate, single event
(my_ev <- sample(events_to_keep, 1))
(my_rep <- sample(100, 1))
quantifs_filtered <- sim_quantifs[[my_rep]]
mat_sf_expression <- sim_mat_sf[[my_rep]]
y <- quantifs_filtered[quantifs_filtered$event_id == my_ev, c("sample_id", "dPSI_nat")] |>
  column_to_rownames("sample_id") |>
  filter(!is.na("dPSI_nat")) |>
  as.matrix()

# get x data
x <- mat_sf_expression[rownames(y),]

n <- nrow(x)
train <- sample(n, round(.7*n))

fit <- cv.glmnet(x = x[train,], y = y[train,], type.measure = "mse",
                 nfolds = 20, alpha = 1, intercept = FALSE)
plot(fit)

predict(fit, newx = x[-train,], s = "lambda.1se") |>
  as.data.frame() |>
  as_tibble(rownames = "sample_id") |>
  rename(predicted = lambda.1se) |>
  add_column(measured = y[-train]) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = measured, y = predicted))





get_rsq_tpr_fdr <- function(.intercept, .unit){
  map(sample(events_to_keep, 10),
      \(my_ev){
        reg_lasso_coefs <- map(sample(100, 10),
                               ~ regression_wrapper(my_ev = my_ev, regression_method = "lasso",
                                                    column = .unit,
                                                    shuffle = FALSE,
                                                    mat_sf_expression = sim_mat_sf[[.x]],
                                                    quants = sim_quantifs[[.x]],
                                                    intercept = .intercept) |>
                                 pluck("coefs_sf") |> 
                                 add_column(sim_rep = paste0("sim_", .x))) |>
          list_rbind()
        
        
        true_coefs <- imap(sim_true_coefs,
                           ~ .x |>
                             filter(event_id == my_ev) |>
                             add_column(sim_rep = paste0("sim_", .y))) |>
          list_rbind()
        
        
        
        # Compare to true known coefficients
        comp_to_true <- left_join(reg_lasso_coefs |>
                                    filter(transcript_id != "(Intercept)"),
                                  true_coefs |>
                                    select(-event_id),
                                  by = c("transcript_id", "sim_rep"))
        
        # Get metrics, 1 Rsquared
        rsquared <- map_dbl(unique(comp_to_true$sim_rep),
                ~ lm(s1 ~ true_coef, data = comp_to_true |>
                       filter(transcript_id != "(Intercept)", sim_rep == .x)) |>
                  summary() |>
                  {\(x) x$adj.r.squared}())
        
        # metrics 2 TPR FDR FPR
        tpr_fdr <- comp_to_true |>
          filter(transcript_id != "(Intercept)") |>
          mutate(s1 = s1 != 0,
                 true_coef = true_coef != 0) |>
          group_by(sim_rep) |>
          summarize(TP = sum(s1 & true_coef),
                    FP = sum(s1 & !true_coef),
                    FN = sum(!s1 & true_coef),
                    TN = sum(!s1 & !true_coef),
                    .groups = 'drop') |>
          mutate(TPR = TP/(TP+FN),
                 FPR = FP/(FP+TN),
                 FDR = FP/(FP+TP))
        
        
        list(rsquared=rsquared, tpr_fdr=tpr_fdr)
        
      })
}


all_sim_fit <- tibble(unit = rep(c("PSI","dPSI_nat","dPSI_logit"), times = 2),
                      intercept = rep(c(TRUE, FALSE), each = 3)) |>
  mutate(res = map2(intercept, unit, get_rsq_tpr_fdr,
                    .progress = TRUE))

res_all_sim_fit <- all_sim_fit |>
  unnest(res) |>
  mutate(rsquared = list_transpose(res)[["rsquared"]],
         tpr_fdr = list_transpose(res)[["tpr_fdr"]]) |>
  select(-res) |>
  unnest(c(rsquared, tpr_fdr))
  



res_all_sim_fit |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = intercept, fill = unit, y = rsquared)) +
  scale_y_continuous(limits = c(0,1))


res_all_sim_fit |>
  select(intercept, unit, TPR, FPR, FDR) |>
  pivot_longer(-c(intercept, unit),
               names_to = "metric") |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = intercept, fill = unit, y = value)) +
  facet_wrap(~metric) +
  xlab(NULL)+ylab(NULL)







# ---- Check permutation test on simulations ----

# Load results ----
# # with lambda min
# perm_sim1 <- qs::qread("data/intermediates/simultation/perm_tests/230203_sim1_regression_permutations_psi.qs")
# perm_sim2 <- qs::qread("data/intermediates/simultation/perm_tests/230203_sim2_regression_permutations_psi.qs")

# with lambda 1se
perm_sim1 <- qs::qread("data/intermediates/simultation/perm_tests/lambda_1se/230203_sim1_regression_permutations_psi_1se.qs")
perm_sim2 <- qs::qread("data/intermediates/simultation/perm_tests/lambda_1se/230203_sim2_regression_permutations_psi_1se.qs")


sim_true_coefs[[1]]


# Plots ----

perm_sim1 |>
  filter(transcript_id != "(Intercept)") |>
  left_join(sim_true_coefs[[1]],
            by = c("event_id", "transcript_id")) |>
  mutate(is_true = true_coef != 0) |> 
  # filter(is_true) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x= coef_effect_size, y = -log10(p_adj),
                 color = is_true, alpha = is_true, size = is_true)) +
  scale_color_manual(values = c("black", "red")) +
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_size_manual(values = c(.5, 2))

# Volcano Plot
perm_sim1 |>
  filter(transcript_id != "(Intercept)") |>
  left_join(sim_true_coefs[[1]],
            by = c("event_id", "transcript_id")) |>
  mutate(is_true = true_coef != 0) |> 
  ggplot(aes(x= coef_effect_size, y = -log10(p_adj),
             color = is_true, alpha = is_true, size = is_true)) +
  theme_classic() +
  geom_point() +
  scale_x_continuous(trans = root_trans(4, signed = TRUE)) +
  xlab(expression(sqrt("effect size", 4))) +
  scale_color_manual(values = c("black", "red")) +
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_size_manual(values = c(.5, 2))

perm_sim2 |>
  filter(transcript_id != "(Intercept)") |>
  left_join(sim_true_coefs[[2]],
            by = c("event_id", "transcript_id")) |>
  mutate(is_true = true_coef != 0) |> 
  ggplot(aes(x= coef_effect_size, y = -log10(p_adj),
             color = is_true, alpha = is_true, size = is_true)) +
  theme_classic() +
  geom_point() +
  scale_x_continuous(trans = root_trans(4, signed = TRUE)) +
  xlab(expression(sqrt("effect size", 4))) +
  scale_color_manual(values = c("black", "red")) +
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_size_manual(values = c(.5, 2))






# 




# ---- For single replicate ----

# Load data ----
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





