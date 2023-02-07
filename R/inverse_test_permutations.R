# same simulations as previous, try the opposite regression

# ie, predict TPM from PSI


library(tidyverse)
library(glmnet)
library(wbData)

gids <- wb_load_gene_ids(281) |>
  add_row(X="0", gene_id = "(Intercept)", symbol ="(Intercept)",
          sequence = "(Intercept)", status="Live",biotype="none",name="(Intercept)")





# ---- Real data ----

# Read data ----

# real
quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs")



sf_tx2g <- sf_expression |>
  select(transcript_id, gene_id, gene_name) |>
  add_row(transcript_id = "(Intercept)",
          gene_id =  "(Intercept)",
          gene_name = "(Intercept)") |>
  distinct()

convert_sf_tx2g <- function(tx_names, warn_missing = TRUE){
  res <- sf_tx2g$gene_id[match(tx_names, sf_tx2g$transcript_id, incomparables = NA)]
  if (warn_missing && any(is.na(res))) {
    warning("i2s: ", sum(is.na(res)), " tx names could not be converted. NA are returned.")
  }
  res
}



events_coordinates <- read_tsv("data/export_for_arman/221111_events_coordinates.tsv")

convert_event2_gene_id <- function(event_ids, warn_missing = TRUE){
  res <- events_coordinates$gene_id[match(event_ids, events_coordinates$event_id, incomparables = NA)]
  if (warn_missing && any(is.na(res))) {
    warning("converts: ", sum(is.na(res)), " event names could not be converted. NA are returned.")
  }
  res
}



# Predict DeltaPSI instead of PSI
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9 for rationale

logit <- function(x){
  stopifnot(all(x >= 0 & x <= 1))
  x[x == 0] <- .01
  x[x == 1] <- .99
  log(x/(1-x))
}

# Prepare quantifs
quantifs_filtered <- quantifs_filtered |>
                      group_by(event_id) |>
                      mutate(dPSI_nat = PSI - mean(PSI),
                             dPSI_logit = logit(PSI) - logit(mean(PSI))) |>
  select(event_id, sample_id, dPSI_nat) |>
  pivot_wider(id_cols = sample_id,
              names_from = "event_id",
              values_from = "dPSI_nat") |>
  column_to_rownames("sample_id")



# Prepare SF
sf_expression <- sf_expression |>
  mutate(logTPM = log(TPM + 1))


all_tx <- unique(sf_expression$transcript_id)



###

# Play with examples ----
(my_tx <- sample(all_tx, 1))

sf_expression |>
  filter(transcript_id == my_tx) |>
  # select(transcript_id, sample_id, TPM, logTPM) |>
  pivot_longer(ends_with("TPM"),
               names_to = "type",
               values_to = "expr") |>
  ggplot() +
  theme_classic() +
  facet_wrap(~type, scales = "free_y") +
  geom_col(aes(x = sample_id, y = expr, fill = sample_id, group = sample_id),
           position = "dodge", show.legend = FALSE, color = 'white') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
  hues::scale_fill_iwanthue()


# regression
y <- sf_expression[sf_expression$transcript_id == my_tx, c("sample_id", "logTPM")] |>
  column_to_rownames("sample_id") |>
  filter(!is.na("logTPM")) |>
  as.matrix()

# get x data
x <- makeX(quantifs_filtered[rownames(y),], na.impute = TRUE)
x <- scale(x, center = TRUE, scale = TRUE)

n <- nrow(x)
train <- sample(n, round(.7*n))

fit <- cv.glmnet(x = x[train,], y = y[train,], type.measure = "mse", nfolds = 20,
                 alpha = 1, scale = FALSE)
plot(fit)




regression_wrapper <- function(my_tx, shuffle, quants, sf_expression){
  
  y <- sf_expression[sf_expression$transcript_id == my_tx, c("sample_id", "logTPM")] |>
    column_to_rownames("sample_id") |>
    filter(!is.na("logTPM")) |>
    as.matrix()
  
  # get x data
  x <- makeX(quants[rownames(y),], na.impute = TRUE)
  x <- scale(x, center = TRUE, scale = TRUE)
  
  n <- nrow(x)
  train <- sample(n, round(.7*n))
  
  fit <- cv.glmnet(x = x[train,], y = y[train,], type.measure = "mse", nfolds = 20,
                   alpha = 1, scale = FALSE)
  
  
  # Estimate on test data
  prediction_on_test <- predict(fit, newx = x[-train,], s = "lambda.1se") |>
    as.data.frame() |>
    as_tibble(rownames = "sample_id") |>
    rename(predicted = lambda.1se) |>
    add_column(measured = y[-train])
  
  rsquare <- summary(lm(predicted ~ measured, data = prediction_on_test))$adj.r.squared
  
  coefs_sf <- coef(fit, s = "lambda.1se") |>
    as.matrix() |>
    as_tibble(rownames = "event_id")
  
  
  
  list(rsquare = rsquare, nb_coefs = sum(coefs_sf$s1 != 0),
       prediction_on_test = prediction_on_test,
       coefs_sf = coefs_sf, fit = fit)
}
sparse_regression <- function(my_tx, shuffle = FALSE, quants, mat_sf_expression){
  possibly(regression_wrapper,
           otherwise = list(NA))(my_tx, shuffle, quants, mat_sf_expression)
}




regression_permutations <- expand_grid(transcript_id = all_tx |> sample(10),
                                       shuffle = c(FALSE, rep(TRUE, 20))) |>
  mutate(res = map2(transcript_id, shuffle,
                    \(.transcript_id, .shuffle){
                      regression_wrapper(.transcript_id, .shuffle, quantifs_filtered, sf_expression)
                      },
                    .progress = TRUE))

regression_permutations <- regression_permutations |>
  mutate(
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)

regression_permutations |>
  ggplot() +
  theme_classic() +
  facet_wrap(~shuffle) +
  geom_point(aes(x = nb_coefs, y = rsquare))




regression_permutations |>
  select(transcript_id, shuffle, rsquare, nb_coefs) |>
  ggplot() +
  theme_classic() +
  geom_jitter(aes(x = nb_coefs, y = rsquare), alpha = .2) +
  facet_wrap(~shuffle)





perm_effect_size <- regression_permutations |>
  select(real_transcript_id = transcript_id, shuffle, coefs_sf) |>
  unnest(coefs_sf) |>
  rename(event_id = transcript_id,
         transcript_id = real_transcript_id) |>
  group_by(shuffle, transcript_id, event_id) |>
  summarize(mean_s = mean(s1),
            .groups = "drop") |>
  pivot_wider(c(transcript_id, event_id),
              names_from = "shuffle",
              values_from = "mean_s") |>
  rename(coef = `FALSE`,
         mean_null = `TRUE`) |>
  mutate(coef_effect_size = coef - mean_null)




perm_p_val <- regression_permutations |>
  select(real_transcript_id = transcript_id, shuffle, coefs_sf) |>
  unnest(coefs_sf) |>
  rename(event_id = transcript_id,
         transcript_id = real_transcript_id) |>
  pivot_wider(c(transcript_id, event_id),
              names_from = "shuffle",
              values_from = "s1",
              values_fn = list) |>
  rename(coef = `FALSE`,
         coefs_null = `TRUE`) |>
  rowwise() |>
  mutate(p_val = sum(abs(coefs_null) >= abs(coef))/length(coefs_null)) |>
  ungroup() |>
  mutate(p_adj = p.adjust(p_val, "BH")) |>
  select(transcript_id, event_id, p_val, p_adj)



perm_res <- left_join(perm_effect_size,
                      perm_p_val,
                      by = c("event_id", "transcript_id")) |>
  mutate(sf_id = convert_sf_tx2g(transcript_id),
         target_id = convert_event2_gene_id(event_id),
         sf_name = i2s(sf_id, gids, warn_missing = TRUE),
         target_name = i2s(target_id, gids, warn_missing = TRUE))

perm_res2 <- perm_res |> 
  mutate(p_adj = p_adj + 1e-5)
perm_res |>
  filter(transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x= coef_effect_size, y = -log10(p_adj),
                 color = p_adj <= 0.1)) +
  scale_x_continuous(trans = root_trans(4, signed = TRUE)) +
  xlab(expression(sqrt("effect size", 4))) +
  # ggrepel::geom_text_repel(aes(x= coef_effect_size, y = -log10(p_adj),
  #                              label = paste0(sf_name,"_",target_name)),
  #                          data = perm_res2 |> filter(p_adj <= 0.1)) +
  scale_color_manual(values = c("black", "red"))





# ~~~ ----

# ---- Simulations ----

# Read data ----

# real
quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs")

# simulation
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

# Prepare quantifs
sim_quantifs <- map(sim_quantifs,
                    ~ .x |>
                      group_by(event_id) |>
                      mutate(dPSI_nat = PSI - mean(PSI),
                             dPSI_logit = logit(PSI) - logit(mean(PSI))) |>
                      select(event_id, sample_id, dPSI_nat) |>
                      pivot_wider(id_cols = sample_id,
                                  names_from = "event_id",
                                  values_from = "dPSI_nat") |>
                      column_to_rownames("sample_id"))


# Prepare SF
sim_sf <- map(sim_sf,
              ~ .x |>
                mutate(logTPM = log(TPM + 1)))


all_tx <- unique(sim_sf[[1]]$transcript_id)



###

# Play with examples ----
(my_tx <- sample(all_tx, 1))

sim_sf[[37]] |>
  filter(transcript_id == my_tx) |>
  # select(transcript_id, sample_id, TPM, logTPM) |>
  pivot_longer(ends_with("TPM"),
               names_to = "type",
               values_to = "expr") |>
  ggplot() +
  theme_classic() +
  facet_wrap(~type, scales = "free_y") +
  geom_col(aes(x = sample_id, y = expr, fill = sample_id, group = sample_id),
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
(my_tx <- sample(all_tx, 1))
(my_rep <- sample(100, 1))
quantifs_filtered <- sim_quantifs[[my_rep]]
sf_expression <- sim_sf[[my_rep]]


# regression
y <- sf_expression[sf_expression$transcript_id == my_tx, c("sample_id", "logTPM")] |>
  column_to_rownames("sample_id") |>
  filter(!is.na("logTPM")) |>
  as.matrix()

# get x data
x <- makeX(quantifs_filtered[rownames(y),], na.impute = TRUE)
x <- scale(x, center = TRUE, scale = TRUE)

n <- nrow(x)
train <- sample(n, round(.7*n))

fit <- cv.glmnet(x = x[train,], y = y[train,], type.measure = "mae",
                 nfolds = 20, alpha = 1, scale = FALSE, intercept = FALSE)
plot(fit)









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
