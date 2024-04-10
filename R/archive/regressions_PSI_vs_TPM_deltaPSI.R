
# Perform regression with different methods. Combination of explore_lasso_on_PSI_vs_TPM and explore_adaptive_lasso_on_PSI_vs_TPM

# Same as others, but for LASSO and adaptive LASSO, predict the deltaPSI instead of the PSI.
# What's more, try both with deltaPSI on natural scale, or taking the logit
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9


# Inits ----

library(tidyverse)
library(glmnet)
library(wbData)


tx2g <- wb_load_tx2gene(281)
gids <- wb_load_gene_ids(281) |>
  add_row(X="0", gene_id = "(Intercept)", symbol ="(Intercept)",
          sequence = "(Intercept)", status="Live",biotype="none",name="(Intercept)")


source("R/regression_functions.R")


# Read data ----

# quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
#   mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) 
quantifs <- qs::qread("data/intermediates/230117_quantifs_filtered.qs")

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

# sf_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") |>
#   filter(gene_id %in% putative_splice_factors) |>
#   mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE))

sf_expression <- qs::qread("data/intermediates/230117_tx_filtered.qs") |>
  filter(gene_id %in% putative_splice_factors) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE))

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



# Filter events ----

# Keep only if enough reads supporting measure


quantifs |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads),color='grey', bins = 200) +
  geom_vline(aes(xintercept = 20), linetype = 'dashed') +
  scale_x_log10()


quantifs_filtered_n_reads <- quantifs |>
  filter(nb_reads > 20)


# filter based on nb of samples that event was measured in
quantifs_filtered_n_reads |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  ggplot() + theme_classic() +
  geom_point(aes(x = nb_neurons, y = nb_samples), alpha = .2) +
  xlab("Number of neurons") + ylab("Number of samples") +
  geom_hline(aes(yintercept = 120.5), color = 'darkred') +
  geom_vline(aes(xintercept = 32.5), color = 'darkred')



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
quantifs_filtered_nsamples |>
  filter(!is.na(PSI)) |>
  group_by(neuron_id, event_id) |>
  summarize(mean_PSI  = mean(PSI)) |>
  group_by(event_id) |>
  summarize(mean_PSI_btw_neurs  = mean(mean_PSI),
            sd_PSI_btw_neurs = sd(mean_PSI),
            .groups = 'drop') |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = mean_PSI_btw_neurs, y = sd_PSI_btw_neurs)) +
  geom_hline(aes(yintercept = 0.05), color = 'grey')

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


quantifs_filtered |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads),color='grey', bins = 100) +
  scale_x_log10()



# # save filtered objects to generate simulations from
# qs::qsave(quantifs_filtered, "data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
# qs::qsave(sf_expression, "data/intermediates/simultation/230206_preprocessed_sf_expression.qs")






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

# Play with examples ----
# SE_136, SE_303 are good. 1060 particularly
(my_ev <- sample(events_to_keep, 1))
# my_ev <- "SE_1046"

events_coordinates |> filter(event_id== my_ev) |> pull(gene_id) |> i2s(gids)

quantifs_filtered |>
  ungroup() |> select(-nb_reads) |>
  filter(event_id == my_ev) |>
  pivot_longer(-ends_with("_id"),
               names_to = "type",
               values_to = "PSI") |>
  ggplot() +
  theme_classic() +
  facet_wrap(~type, scales = "free_y") +
  geom_col(aes(x = neuron_id, y = PSI, fill = neuron_id, group = sample_id),
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


# try replicating
out <- list()
for(i in 1:5){
  reg_lasso <- regression_lasso_psi(my_ev, mat_sf_expression, quantifs_filtered)
  reg_adalasso <- regression_adalasso_psi(my_ev, mat_sf_expression, quantifs_filtered)
  
  reg_lasso_nat <- regression_lasso_dpsi_nat(my_ev, mat_sf_expression, quantifs_filtered)
  reg_adalasso_nat <- regression_adalasso_dpsi_nat(my_ev, mat_sf_expression, quantifs_filtered)
  
  reg_lasso_logit <- regression_lasso_dpsi_logit(my_ev, mat_sf_expression, quantifs_filtered)
  reg_adalasso_logit <- regression_adalasso_dpsi_logit(my_ev, mat_sf_expression, quantifs_filtered)
  
  
  out[[i]] <- tibble(fit = list(reg_lasso, reg_lasso_nat, reg_lasso_logit,
                    reg_adalasso, reg_adalasso_nat, reg_adalasso_logit),
         method = c("reg_lasso", "reg_lasso_nat", "reg_lasso_logit",
                    "reg_adalasso", "reg_adalasso_nat", "reg_adalasso_logit"),
         rsquare = map_dbl(fit, ~ pluck(.x, "rsquare")))
}
out |>
  list_rbind() |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = method, y = rsquare))





# coefficients
reg_lasso$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  left_join(sf_tx2g,
            by = "transcript_id") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = gene_name, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


reg_adalasso$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  left_join(sf_tx2g,
            by = "transcript_id") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = gene_name, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 

reg_lasso_nat$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  left_join(sf_tx2g,
            by = "transcript_id") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = gene_name, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


reg_adalasso_nat$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  left_join(sf_tx2g,
            by = "transcript_id") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = gene_name, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 

reg_lasso_logit$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  left_join(sf_tx2g,
            by = "transcript_id") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = gene_name, y = s1)) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


reg_adalasso_logit$coefs_sf |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  left_join(sf_tx2g,
            by = "transcript_id") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = gene_name, y = s1)) +
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




# First pass ----
# make a first round of regressions

first_pass <- expand_grid(event_id = unique(quantifs_filtered$event_id),
                              method = c("lasso", "adalasso"),
                              column = c("PSI", "dPSI_nat", "dPSI_logit")) |>
  # slice_sample(n = 2) |>
  mutate(res = pmap(list(event_id, method, column),
                    \(event_id, method, column) sparse_regression(event_id, method, column,
                                                                  shuffle = FALSE, mat_sf_expression, quantifs_filtered),
                    .progress = TRUE),
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)


####  Using permutation testing ----


regression_permutations <- expand_grid(event_id = sample(unique(quantifs_filtered$event_id), 20),
                              method = c("lasso"),
                              column = c("dPSI_nat"),
                              shuffle = c(FALSE, rep(TRUE, 10))) |>
  mutate(res = pmap(list(event_id, method, column, shuffle),
                    \(event_id, method, column, shuffle) sparse_regression(event_id, method, column,
                                                                  shuffle, mat_sf_expression, quantifs_filtered),
                    .progress = TRUE),
         fit = map(res, ~pluck(.x, "fit", .default = NA)),
         rsquare = map_dbl(res, ~pluck(.x, "rsquare", .default = NA_real_)),
         nb_coefs = map_int(res, ~pluck(.x, "nb_coefs", .default = NA_integer_)),
         prediction_on_test = map(res, ~pluck(.x, "prediction_on_test", .default = tibble())),
         coefs_sf = map(res, ~pluck(.x, "coefs_sf", .default = tibble()))) |>
  select(-res)

# qs::qsave(regression_permutations, "data/intermediates/230203_regression_permutations_test_dPSI.qs")

regression_permutations <- qs::qread("data/intermediates/230202_regression_permutations.qs")

regression_permutations |>
  select(event_id, shuffle, rsquare, nb_coefs) |>
  ggplot() +
  theme_classic() +
  geom_jitter(aes(x = nb_coefs, y = rsquare), alpha = .2) +
  facet_wrap(~shuffle)





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
          by = c("event_id", "transcript_id")) |>
  mutate(sf_id = convert_sf_tx2g(transcript_id),
         target_id = convert_event2_gene_id(event_id),
         sf_name = i2s(sf_id, gids, warn_missing = TRUE),
         target_name = i2s(target_id, gids, warn_missing = TRUE))

# qs::qsave(perm_res, "data/intermediates/230202_permutations_results_from_ruddle.qs")

perm_res_rud <- qs::qread("data/intermediates/230202_permutations_results_from_ruddle.qs")

perm_res |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x= coef_effect_size, y = -log10(p_adj)), alpha = .2) +
  # scale_x_continuous(limits = c(-.05,.05)) +
  NULL

# On log scale
perm_res |>
  ggplot(aes(x= log10(abs(coef_effect_size)), y = -log10(p_adj))) +
  theme_classic() +
  geom_point(aes(color = as.factor(sign(coef_effect_size))))
  # ggrepel::geom_text_repel(aes(label = paste0(sf_name,"_",target_name)),
  #                          data = perm_res |> filter(p_adj < 1))

root_breaks <- function(n = 10, exponent, signed){
  n_default <- n
  function(x, n = n_default){
    if(signed){
      min <- min(abs(x))
      max <- max(abs(x))
    } else{
      min <- min(x)
      max <- max(x)
    }
    
    by <- (max - min)/n
    breaks <- seq(min, 1.1*max^(1/exponent), by = by)^exponent
    breaks <- round(breaks, digits = 2)
    
    if(signed){
      return(c(-breaks, 0, breaks))
    } else{
      return(c(0,breaks))
    }
  }
}

root_trans <- function(exponent = 2, signed = FALSE){
  if(signed){
    scales::trans_new("root_trans",
                      \(x) sign(x)*abs(x)^(1/exponent),
                      \(x) sign(x)*abs(x)^exponent,
                      breaks = root_breaks(exponent = exponent, signed = signed),
                      domain = c(-Inf,Inf))
  } else{
    scales::trans_new("root_trans",
                      \(x) x^(1/exponent),
                      \(x) x^exponent,
                      breaks = root_breaks(exponent = exponent, signed = signed),
                      domain = c(0,Inf))
  }
}




perm_res |>
  filter(p_adj < 1) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = coef_effect_size), bins = 100, color = 'white') +
  scale_x_continuous(trans = root_trans(4, signed = TRUE))


# Volcano Plot
perm_res |>#slice_sample(n=100) |> 
  mutate(selected = abs(coef_effect_size) >= 0.005 & p_adj <= 0.1) |>
  ggplot(aes(x= coef_effect_size, y = -log10(p_adj), color = selected)) +
  theme_classic() +
  geom_point(alpha = .2) +
  scale_x_continuous(trans = root_trans(4, signed = TRUE)) +
  xlab(expression(sqrt("effect size", 4))) +
  theme(axis.text = element_text(size = 7)) +
  geom_vline(aes(xintercept = .005)) +
  geom_vline(aes(xintercept = -.005)) +
  geom_hline(aes(yintercept = -log10(.1)))

perm_res_rud |>#slice_sample(n=10000) |> 
  mutate(selected = abs(coef_effect_size) >= 0.005 & p_adj <= 0.1) |>
  ggplot(aes(x= coef_effect_size, y = -log10(p_adj), color = selected)) +
  theme_classic() +
  geom_point(alpha = .2) +
  scale_x_continuous(trans = root_trans(4, signed = TRUE)) +
  xlab(expression(sqrt("effect size", 4))) +
  theme(axis.text = element_text(size = 7)) +
  geom_vline(aes(xintercept = .005)) +
  geom_vline(aes(xintercept = -.005)) +
  geom_hline(aes(yintercept = -log10(.1)))


interactions_rud <- perm_res_rud |>
  mutate(selected = abs(coef_effect_size) >= 0.005 & p_adj <= 0.1) |>
  filter(selected) |>
  mutate(interaction = paste0(sf_name,"_",target_name)) |>
  pull(interaction)


interactions_dpsi <- perm_res |>
  mutate(selected = abs(coef_effect_size) >= 0.005 & p_adj <= 0.1) |>
  filter(selected) |>
  mutate(interaction = paste0(sf_name,"_",target_name)) |>
  pull(interaction)


eulerr::euler(list(dpsi = unique(interactions_dpsi),
                   rud = unique(interactions_rud))) |>
  plot(quantities = TRUE)


########
sparse_regression(event_id, method, column,
                  shuffle = FALSE, mat_sf_expression, quantifs_filtered)

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









# Old code ----

first_pass_regression <- quantifs_filtered |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  arrange(desc(nb_samples), desc(nb_neurons)) |>
  # slice_head(n = 2) |>
  mutate(fit_lasso = map(event_id,
                         \(ev) possibly(regression_lasso_psi,
                                        otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                         .progress = TRUE),
         fit_adalasso = map(event_id,
                            \(ev) possibly(regression_adalasso_psi,
                                           otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                            .progress = TRUE),
         fit_lasso_nat = map(event_id,
                         \(ev) possibly(regression_lasso_dpsi_nat,
                                        otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                         .progress = TRUE),
         fit_adalasso_nat = map(event_id,
                            \(ev) possibly(regression_adalasso_dpsi_nat,
                                           otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                            .progress = TRUE),
         fit_lasso_logit = map(event_id,
                         \(ev) possibly(regression_lasso_dpsi_logit,
                                        otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                         .progress = TRUE),
         fit_adalasso_logit = map(event_id,
                            \(ev) possibly(regression_adalasso_dpsi_logit,
                                           otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                            .progress = TRUE))
# qs::qsave(first_pass_regression, "data/intermediates/230201_fits_for_quantifs_dpsi_cache.qs")


# first_pass_regression <- qs::qread("data/intermediates/230125_fits_for_quantifs_4methods_cache.qs")


# Compare results on first pass ----

patchwork::wrap_plots(
  ggplot(first_pass_regression |>
           mutate(rsquare = map_dbl(fit_lasso, pluck("rsquare"))),
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(rsquare = map_dbl(fit_adalasso, pluck("rsquare"))),
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(rsquare = map_dbl(fit_lasso_nat, pluck("rsquare"))),
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(rsquare = map_dbl(fit_adalasso_nat, pluck("rsquare"))),
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(rsquare = map_dbl(fit_lasso_logit, pluck("rsquare"))),
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(rsquare = map_dbl(fit_adalasso_logit, pluck("rsquare"))),
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ncol = 2
)

patchwork::wrap_plots(
  ggplot(first_pass_regression |>
           mutate(nb_coefs = map_int(fit_lasso, ~ pluck(., "nb_coefs"))),
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(nb_coefs = map_int(fit_adalasso, ~ pluck(., "nb_coefs"))),
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(nb_coefs = map_int(fit_lasso_nat, ~ pluck(., "nb_coefs"))),
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(nb_coefs = map_int(fit_adalasso_nat, ~ pluck(., "nb_coefs"))),
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(nb_coefs = map_int(fit_lasso_logit, ~ pluck(., "nb_coefs"))),
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(nb_coefs = map_int(fit_adalasso_logit, ~ pluck(., "nb_coefs"))),
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ncol = 2
)




first_pass_regression |>
  select(event_id, ends_with("_rsquare")) |>
  pivot_longer(-event_id, names_to = "method", values_to = "rsquare", names_pattern = "(a?d?a?lasso)_rsquare") |>
  ggplot(aes(x = method, y = rsquare)) +
  theme_classic() +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm()

first_pass_regression |>
  select(event_id, ends_with("_coefs_sf")) |>
  pivot_longer(-event_id, names_to = "method", values_to = "coefs", names_pattern = "(a?d?a?lasso)_coefs_sf") |>
  unnest(coefs) |>
  filter(s1 != 0,
         transcript_id != "(Intercept)") |>
  ggplot(aes(x = transcript_id, y = s1)) +
  theme_classic() +
  geom_boxplot() +
  facet_wrap(~ method)
ggbeeswarm::geom_beeswarm()


full_join(first_pass_regression |>
            select(event_id, lasso_coefs_sf) |>
            unnest(lasso_coefs_sf) |>
            rename(lasso_coef = s1),
          first_pass_regression |>
            select(event_id, adalasso_coefs_sf) |>
            unnest(adalasso_coefs_sf) |>
            rename(adalasso_coef = s1)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = lasso_coef, y = adalasso_coef))






# Check known regulated events ----

sf_targets <- jsonlite::read_json("data/export_for_arman/sf_targets_v2.json")

sf2target <- tibble(sf_id = map_chr(sf_targets, \(x) x[["SF"]]),
                    target_id = map(sf_targets, \(x) x[["targets"]])) |>
  mutate(sf_id = if_else(sf_id == "mec-8 ad", "mec-8", sf_id),
         sf_id = map_chr(sf_id,
                         \(x) `if`(startsWith(x, "WBGene"), x, s2i(x, gids, warn_missing = TRUE)))) |>
  unnest_longer(target_id) |>
  filter(target_id != "0") |>
  mutate(sf_name = i2s(sf_id, gids, warn_missing = FALSE),
         target_name = i2s(target_id, gids, warn_missing = TRUE))

# how many of the events have known regulators?
table(events_coordinates$gene_id %in% sf2target$target_id)

events_with_known_regulators <- events_coordinates$event_id[events_coordinates$gene_id %in% sf2target$target_id]

table(first_pass_regression$event_id %in% events_with_known_regulators)

first_pass_filt <- first_pass_regression |>
  filter(event_id %in% events_with_known_regulators,
         event_id %in% events_to_keep)





overlaps_sf <- first_pass_filt |>
  mutate(adalasso_computed_sf = map(adalasso_coefs_sf,
                                    \(x) convert_sf_tx2g(x$transcript_id[x$s1 !=0 & x$transcript_id != "(Intercept)"])),
         lasso_computed_sf = map(lasso_coefs_sf,
                                 \(x) convert_sf_tx2g(x$transcript_id[x$s1 !=0 & x$transcript_id != "(Intercept)"]))) |>
  left_join(events_coordinates |>
              select(event_id, gene_id),
            by = "event_id") |>
  left_join(sf2target,
            by = c(gene_id = "target_id")) |>
  select(-sf_name) |>
  group_by(across(-sf_id)) |>
  summarize(known_sf = list(sf_id),
            .groups = "drop") |>
  mutate(nb_lasso_overlapping_sf = map2_int(lasso_computed_sf, known_sf,
                                            \(x,y) length(intersect(x,y))),
         nb_adalasso_overlapping_sf = map2_int(adalasso_computed_sf, known_sf,
                                               \(x,y) length(intersect(x,y))),
         nb_known_sf = map_int(known_sf, length),
         nb_lasso_computed_sf = map_int(lasso_computed_sf, length),
         nb_adalasso_computed_sf = map_int(adalasso_computed_sf, length))

overlaps_sf |>
  mutate(prop_lasso_overlapping_sf = nb_lasso_overlapping_sf/nb_known_sf,
         prop_adalasso_overlapping_sf = nb_adalasso_overlapping_sf/nb_known_sf) |>
  select(nb_samples, starts_with("prop_")) |>
  pivot_longer(-nb_samples,
               names_to = "method", names_pattern = "prop_(a?d?a?lasso)_overlapping_sf",
               values_to = "prop_overlapping_sf") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nb_samples, y = prop_overlapping_sf), alpha = .2) +
  facet_wrap(~method)

sum(overlaps_sf$nb_known_sf)
sum(overlaps_sf$nb_lasso_computed_sf)
sum(overlaps_sf$nb_lasso_overlapping_sf)


sum(overlaps_sf$nb_known_sf)
sum(overlaps_sf$nb_adalasso_computed_sf)
sum(overlaps_sf$nb_adalasso_overlapping_sf)




# randomize and test
randomized_overlaps <- pbapply::pbreplicate(200,{
  first_pass_filt_rand <- first_pass_filt
  first_pass_filt_rand$coefs_sf <- sample(first_pass_filt$coefs_sf)
  
  xx_rand <- first_pass_filt_rand |>
    mutate(computed_sf = map(coefs_sf,
                             \(x) convert_sf_tx2g(x$transcript_id[x$s1 !=0 & x$transcript_id != "(Intercept)"])
    )) |>
    left_join(events_coordinates |>
                select(event_id, gene_id),
              by = "event_id") |>
    left_join(sf2target,
              by = c(gene_id = "target_id")) |>
    select(-sf_name) |>
    group_by(across(-sf_id)) |>
    summarize(known_sf = list(sf_id),
              .groups = "drop") |>
    mutate(nb_overlapping_sf = map2_int(computed_sf, known_sf,
                                        \(x,y) length(intersect(x,y))),
           nb_known_sf = map_int(known_sf, length),
           nb_computed_sf = map_int(computed_sf, length))
  sum(xx_rand$nb_overlapping_sf)
})

hist(randomized_overlaps, breaks = 50) ; abline(v = sum(overlaps_sf$nb_overlapping_sf), col = 'darkred')

table(randomized_overlaps >= sum(overlaps_sf$nb_overlapping_sf))
5/200
24/200
15/200
0/200



# Check reproducibility ----

# first_pass_regression |>
#   ggplot() + theme_classic() +
#   geom_point(aes(x = nb_neurons, y = nb_samples), alpha = .2) +
#   xlab("Number of neurons") + ylab("Number of samples") +
#   geom_hline(aes(yintercept = 178), color = 'darkred') +
#   geom_vline(aes(xintercept = 45.5), color = 'darkred')
# 
# table(first_pass_regression$nb_samples > 178 &
#         first_pass_regression$nb_neurons > 45)
# 
# events_to_keep2 <- quantifs |>
#   filter(! is.na(PSI)) |>
#   group_by(event_id) |>
#   summarize(nb_samples = n(),
#             nb_neurons = n_distinct(neuron_id)) |>
#   filter(nb_samples > 178,
#          nb_neurons > 45) |>
#   pull(event_id)
# 
# 
# table(unique(quantifs$event_id) %in% events_to_keep2)
# 
# quantifs_filtered2 <- quantifs |>
#   filter(! is.na(PSI)) |>
#   filter(event_id %in% events_to_keep2)



# Note: takes 3h
# replicated_regression <- pbapply::pbreplicate(n = 50,
#                                    expr = map(events_to_keep, possibly(sparse_regression,
#                                                                         otherwise = list(rsquare=NA_real_,nb_coefs=NA_real_))),
#                                    simplify = FALSE)
# qs::qsave(replicated_regression, "data/intermediates/230117_replicated_regression_log_cache_filt.qs")


replicated_regression <- qs::qread("data/intermediates/230117_replicated_regression_log_cache_filt.qs")

replicated_rsquare <- map(replicated_regression,
                          \(replicate) {
                            map_dbl(replicate, \(rep) rep[["rsquare"]]) |>
                              set_names(events_to_keep)
                          }) |>
  set_names(paste0("replicate_", 1:50)) |>
  as_tibble() |>
  add_column(event_id = events_to_keep, .before = 1) |>
  pivot_longer(-event_id,
               names_to = "replicate",
               values_to = "Rsquare_adjusted")

ggplot(replicated_rsquare,
       aes(x = event_id, y = Rsquare_adjusted)) +
  theme_classic() +
  # geom_violin(fill = 'grey95') +
  geom_boxplot(fill = 'grey90') +
  # geom_point() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


replicated_rsquare |>
  group_by(event_id) |>
  summarize(mean = mean(Rsquare_adjusted),
            sd = sd(Rsquare_adjusted)) |>
  ggplot(aes(x = mean, y = sd)) +
  theme_classic() +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1), color = 'grey', linetype = 'dashed')

# recheck 1 value
(lm(measured ~ predicted, data = replicated_regression[[2]][[3]]$prediction_on_test) |>
    summary())$adj.r.squared

replicated_rsquare |>
  filter(replicate == "replicate_2", event_id == events_to_keep[[3]]) |>
  pull(Rsquare_adjusted)


# also look at SF computed
replicated_regression[[1]][[147]][["coefs_sf"]]
replicated_coefs <- imap(replicated_regression,
                         \(replicate, ind) {
                           map2(replicate, events_to_keep,
                                \(rep, .event_id) {
                                  pluck(rep, "coefs_sf", .default = tibble()) |>
                                    add_column(event_id = .event_id)
                                }) |>
                             list_rbind() |>
                             add_column(replicate = paste0("replicate_", ind))
                         }) |>
  list_rbind()


(my_ev <- sample(unique(replicated_coefs$event_id), 1))
# xx <- replicated_coefs |>
#   filter(event_id == my_ev,
#          s1 != 0) |>
#   left_join(sf_tx2g,
#             by = "transcript_id") |>
#   ggplot() +
#   theme_classic() +
#   geom_hline(aes(yintercept = 0), color = 'grey') +
#   geom_point(aes(x = gene_name, y = s1, color = replicate), show.legend = FALSE) +
#   theme(
#     axis.text.x = element_text(
#       angle = 90,
#       hjust = 1,
#       vjust = 0.5
#     ))
# 
# plotly::ggplotly(xx)



# Count intersection size
intersected_coefs_single_event <- replicated_coefs |>
  filter(event_id == my_ev,
         s1 != 0) |>
  group_by(transcript_id) |>
  summarize(nb_times_non_zero = n(),
            .groups = "drop") |>
  left_join(sf_tx2g,
            by = "transcript_id")

# find known sf for that event
my_target <- events_coordinates |>
  filter(event_id == my_ev) |> 
  pull(gene_id)
my_sf_ids <- sf2target$sf_name[sf2target$target_id == my_target]

intersected_coefs_single_event |>
  group_by(gene_id) |>
  mutate(mx = max(nb_times_non_zero),
         sm = sum(nb_times_non_zero)) |>
  ungroup() |>
  arrange(desc(mx), desc(sm)) |>
  select(-mx, -sm) |>
  mutate(gene_name = fct_inorder(gene_name),
         known_sf = gene_name %in% my_sf_ids) |>
  # slice_head(prop = .5) |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = gene_name, y = nb_times_non_zero, color = known_sf)) +
  scale_color_manual(values = c("black", "red")) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
  geom_hline(aes(yintercept = 20), color = 'darkred', linetype = "dotted")



# Compare with known ----

intersected_coefs <- replicated_coefs |>
  group_by(transcript_id, event_id) |>
  summarize(nb_intersections = sum(s1 != 0),
            nb_tests = n(),
            .groups = "drop") |>
  # robustness criterion
  filter(nb_intersections >= 20) |>
  # add event info
  left_join(events_coordinates |>
              select(event_id, gene_id),
            by = "event_id") |>
  rename(target_id = gene_id,
         computed_sf_tx_id = transcript_id) |>
  # nest computed sf list
  mutate(computed_sf_gene_id = convert_sf_tx2g(computed_sf_tx_id)) |>
  group_by(event_id, target_id) |>
  nest(computed = c(computed_sf_tx_id, computed_sf_gene_id, nb_intersections, nb_tests)) |>
  ungroup() |>
  # add and nest known sf list
  left_join(sf2target,
            by = "target_id") |>
  group_by(event_id, target_id) |>
  nest(known = c(sf_id, sf_name))

# find intersection size of computed and known SF
overlap_intersected_coefs <- intersected_coefs |>
  select(event_id, target_id, target_name) |>
  add_column(nb_sf_overlap = map2_int(intersected_coefs$computed, intersected_coefs$known,
                                      ~ length(intersect(.x$computed_sf_gene_id, .y$sf_id))),
             nb_sf_computed = map_int(intersected_coefs$computed, ~length(.x$computed_sf_gene_id) - 1),
             nb_sf_known = map_int(intersected_coefs$known, ~length(.x$sf_id)))

sum(overlap_intersected_coefs$nb_sf_computed)
sum(overlap_intersected_coefs$nb_sf_overlap)
sum(overlap_intersected_coefs$nb_sf_known)
25/482
#> 5% no filtering
178/1768
#> 10% after filtering


rep_overlap <- pbapply::pbreplicate(500,
                                    sum(map2_int(sample(intersected_coefs$computed), sample(intersected_coefs$known),
                                                 ~ length(intersect(.x$computed_sf_gene_id, .y$sf_id))))
)

hist(rep_overlap, breaks = 30,
     xlab = "Number of overlapping interactions under randomization",
     main = NULL); abline(v = 178, col = 'red')
table(rep_overlap >= sum(overlap_intersected_coefs$nb_sf_overlap))





# Select only good fits ----
hist(replicated_rsquare$Rsquare_adjusted, breaks = 60); abline(v = .51, col = 'red')
hist(replicated_rsquare$Rsquare_adjusted, breaks = 800, xlim = c(.48,.52)); abline(v = .51, col = 'red')
table(replicated_rsquare$Rsquare_adjusted > .51)

pluck(replicated_regression[[45]][[25]],"coefs_sf", .default = tibble()) |>
  add_column(event_id = "SE111")
length(replicated_regression[[45]])
events_to_keep[24:26]



filt_replicated_coefs <- replicated_coefs |>
  left_join(replicated_rsquare,
            by = c("event_id", "replicate")) |>
  filter(Rsquare_adjusted > 0.51)


#~ Compare with known ----

filt_intersected_coefs <- filt_replicated_coefs |>
  group_by(transcript_id, event_id) |>
  summarize(nb_intersections = sum(s1 != 0),
            nb_tests = n(),
            .groups = "drop") |>
  # robustness criterion (20/50, i.e. 40% of tests)
  filter(nb_intersections/nb_tests >= .4) |>
  # add event info
  left_join(events_coordinates |>
              select(event_id, gene_id),
            by = "event_id") |>
  rename(target_id = gene_id,
         computed_sf_tx_id = transcript_id) |>
  # nest computed sf list
  mutate(computed_sf_gene_id = convert_sf_tx2g(computed_sf_tx_id)) |>
  group_by(event_id, target_id) |>
  nest(computed = c(computed_sf_tx_id, computed_sf_gene_id, nb_intersections, nb_tests)) |>
  ungroup() |>
  # add and nest known sf list
  left_join(sf2target,
            by = "target_id") |>
  group_by(event_id, target_id) |>
  nest(known = c(sf_id, sf_name))

# find intersection size of computed and known SF
filt_overlap_intersected_coefs <- filt_intersected_coefs |>
  select(event_id, target_id, target_name) |>
  add_column(nb_sf_overlap = map2_int(filt_intersected_coefs$computed, filt_intersected_coefs$known,
                                      ~ length(intersect(.x$computed_sf_gene_id, .y$sf_id))),
             nb_sf_computed = map_int(filt_intersected_coefs$computed, ~length(.x$computed_sf_gene_id) - 1),
             nb_sf_known = map_int(filt_intersected_coefs$known, ~length(.x$sf_id)))

sum(filt_overlap_intersected_coefs$nb_sf_overlap)
sum(filt_overlap_intersected_coefs$nb_sf_known)
10/81
#> 12%
50/376
#> 13% after filtering

rep_overlap <- replicate(500,
                         sum(map2_int(sample(filt_intersected_coefs$computed), sample(filt_intersected_coefs$known),
                                      ~ length(intersect(.x$computed_sf_gene_id, .y$sf_id))))
)

hist(rep_overlap, breaks = 10,
     xlab = "Number of overlapping interactions under randomization",
     main = NULL); abline(v = 50, col = 'red')

#> not better





# Find stable coefficients ----
coefs_stability <- replicated_coefs |>
  group_by(event_id, transcript_id) |>
  summarize(median = median(s1),
            mean = mean(s1),
            mad = mad(s1),
            sd = sd(s1),
            .groups = "drop")

ggplot(coefs_stability) +
  theme_classic() +
  geom_point(aes(x = median, y = mad))

ggplot(coefs_stability) +
  theme_classic() +
  geom_point(aes(x = mean, y = sd))

plot(x = coefs_stability$mean, y = coefs_stability$sd)
plot(x = coefs_stability$median, y = coefs_stability$mad)








