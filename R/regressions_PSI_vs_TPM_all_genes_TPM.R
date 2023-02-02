
# Perform regression with different methods. Combination of explore_lasso_on_PSI_vs_TPM and explore_adaptive_lasso_on_PSI_vs_TPM
# this time, use all genes, see if we extract known SFs


# Inits ----

library(tidyverse)
library(glmnet)
library(wbData)


tx2g <- wb_load_tx2gene(281)
gids <- wb_load_gene_ids(281)

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
  # filter(gene_id %in% putative_splice_factors) |>
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


# Filter events ----

# filter based on nb of samples that event was measured in
quantifs |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  ggplot() + theme_classic() +
  geom_point(aes(x = nb_neurons, y = nb_samples), alpha = .2) +
  xlab("Number of neurons") + ylab("Number of samples") +
  geom_hline(aes(yintercept = 120.5), color = 'darkred') +
  geom_vline(aes(xintercept = 32.5), color = 'darkred')



events_to_keep_n_samples <- quantifs |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  filter(nb_samples > 120,
         nb_neurons > 32) |>
  pull(event_id)


quantifs_filtered_nsamples <- quantifs |>
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

quantifs_filtered_nsamples |>
  filter(!is.na(PSI)) |>
  group_by(neuron_id, event_id) |>
  summarize(mean_PSI  = mean(PSI),
            sd_PSI = sd(PSI)) |>
  group_by(event_id) |>
  summarize(mean_PSI_btw_neurs  = mean(mean_PSI),
            sd_PSI_btw_neurs = sd(mean_PSI),
            .groups = 'drop') |>
  arrange(sd_PSI_btw_neurs)

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

quantifs_filtered <- quantifs |>
  filter(event_id %in% events_to_keep) |>
  filter(! is.na(PSI))




# Make binary or ternary events ----
quantifs_filtered |>
  ggplot() +
  theme_classic() +
  geom_freqpoly(aes(x = PSI), bins = 100 ) +
  geom_vline(aes(xintercept = .1)) +
  geom_vline(aes(xintercept = .9)) +
  geom_vline(aes(xintercept = .25)) +
  geom_vline(aes(xintercept = .75))

quantifs_filtered$inclusion <- cut(quantifs_filtered$PSI,
                                   breaks = c(0, .25, .75, 1),
                                   labels = c("excl","mix","incl"),
                                   include.lowest = TRUE)
quantifs_filtered_binary <- quantifs_filtered |>
  filter(inclusion != "mix")


# Filter tx expression ----

# xx <- sf_expression |>
#   group_by(transcript_id) |>
#   summarize(nb_samples = n(),
#             mean = mean(TPM),
#             sd = sd(TPM))
# 
# table(xx$nb_samples, useNA = 'ifany')
# sum(is.na(sf_expression$TPM))
# 
# ggplot(xx) +
#   theme_classic() +
#   geom_point(aes(x = mean, y = sd)) +
#   scale_x_log10() + scale_y_log10()

#> no criterion to filter anything out



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


#~ Refilter ----

# binary
quantifs_filtered_binary |>
  group_by(event_id) |>
  summarize(prop_incl = mean(inclusion == "incl")) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = prop_incl))

quantifs_filtered_binary |>
  group_by(event_id) |>
  summarize(nb_incl = sum(inclusion == "incl"),
            nb_excl = sum(inclusion == "excl")) |>
  pivot_longer(-event_id,
               names_to = "type",
               names_prefix = "nb_",
               values_to = "nb") |>
  ggplot() +
  theme_classic() +
  geom_freqpoly(aes(x = nb, color = type))

quantifs_filtered_binary |>
  group_by(event_id) |>
  summarize(prop_incl = mean(inclusion == "incl"),
            nb_incl = sum(inclusion == "incl"),
            nb_excl = sum(inclusion == "excl")) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nb_incl, y = nb_excl, color = prop_incl))

quantifs_filtered_binary |>
  group_by(event_id) |>
  summarize(prop_incl = mean(inclusion == "incl"),
            nb_incl = sum(inclusion == "incl"),
            nb_excl = sum(inclusion == "excl")) |>
  mutate(bad = prop_incl < .1 | prop_incl > .9) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nb_incl, y = nb_excl, color = bad))

events_to_keep_binary_var <- quantifs_filtered_binary |>
  group_by(event_id) |>
  summarize(prop_incl = mean(inclusion == "incl"),
            nb_incl = sum(inclusion == "incl"),
            nb_excl = sum(inclusion == "excl")) |>
  mutate(bad = prop_incl < .1 | prop_incl > .9) |>
  filter(!bad) |>
  pull(event_id)

quantifs_filtered_binary <- quantifs_filtered_binary |>
  filter(event_id %in% events_to_keep_binary_var)


# same with ternary

quantifs_filtered |>
  group_by(event_id) |>
  summarize(prop_incl = mean(inclusion == "incl"),
            prop_excl = mean(inclusion == "excl"),
            prop_mix = mean(inclusion == "mix")) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = prop_incl, y = prop_excl, color = prop_mix))

get_which_max <- function(x){
  apply(x, 1, \(.row) colnames(x)[which.max(.row)])
}

xx <- quantifs_filtered |>
  group_by(event_id) |>
  summarize(prop_incl = mean(inclusion == "incl"),
            prop_excl = mean(inclusion == "excl"),
            prop_mix = mean(inclusion == "mix"))
xx |>
  mutate(which_max = get_which_max(select(xx, starts_with("prop")))) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = prop_incl, y = prop_excl, color = prop_mix, shape = which_max))

xx |>
  mutate(max = apply(select(xx, starts_with("prop")), 1, max),
         bad = max >= .9) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = prop_incl, y = prop_excl, color = max, shape = bad))

events_to_keep_ternary_var <- xx |>
  mutate(max = apply(select(xx, starts_with("prop")), 1, max),
         bad = max >= .9) |>
  filter(!bad) |>
  pull(event_id)



quantifs_filtered_ternary <- quantifs_filtered |>
  filter(event_id %in% events_to_keep_ternary_var)




# Final list

events_to_keep_all <- events_to_keep |>
  intersect(events_to_keep_binary_var) |>
  intersect(events_to_keep_ternary_var)



# Regression Functions ----

# From https://www.rpubs.com/kaz_yos/alasso
adaLasso <- function(x_cont, y_cont, gamma = 1, type.measure = "mse", nfolds = 10){
  ridge1_cv <- cv.glmnet(x = x_cont, y = y_cont, type.measure = type.measure, nfolds = nfolds, alpha = 0)
  best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
  
  cv.glmnet(x = x_cont, y = y_cont, type.measure = type.measure, nfolds = nfolds, alpha = 1,
            penalty.factor = 1 / (abs(best_ridge_coef)^gamma))
}


regression_lasso <- function(my_ev, mat_sf_expression, quants){
  
  # get y data
  y <- quants[quants$event_id == my_ev, c("sample_id", "PSI")] |>
    column_to_rownames("sample_id") |>
    filter(!is.na(PSI)) |>
    as.matrix()
  
  x <- mat_sf_expression[rownames(y),]
  
  
  n <- nrow(x)
  train <- sample(n, round(.7*n))
  
  
  
  cvfit <- cv.glmnet(x[train,], y[train], nfolds = 20)
  
  # Estimate on test data
  prediction_on_test <- predict(cvfit, newx = x[-train,], s = "lambda.min") |>
    as.data.frame() |>
    as_tibble(rownames = "sample_id") |>
    rename(predicted = lambda.min) |>
    add_column(measured = y[-train])
  
  rsquare <- summary(lm(predicted ~ measured, data = prediction_on_test))$adj.r.squared
  
  coefs_sf <- coef(cvfit, s = "lambda.min") |>
    as.matrix() |>
    as_tibble(rownames = "transcript_id")
  
  nb_coefs <- nrow(coefs_sf |>
                     filter(s1 != 0))
  
  list(rsquare = rsquare, nb_coefs = nb_coefs,
       prediction_on_test = prediction_on_test,
       coefs_sf = coefs_sf)
}


regression_adalasso <- function(my_ev, mat_sf_expression, quants){
  
  # get y data
  y <- quants[quants$event_id == my_ev, c("sample_id", "PSI")] |>
    column_to_rownames("sample_id") |>
    filter(!is.na(PSI)) |>
    as.matrix()
  
  x <- mat_sf_expression[rownames(y),]
  
  n <- nrow(x)
  train <- sample(n, round(.7*n))
  
  
  
  cvfit <- adaLasso(x[train,], y[train], nfolds = 20)
  
  # Estimate on test data
  prediction_on_test <- predict(cvfit, newx = x[-train,], s = "lambda.min") |>
    as.data.frame() |>
    as_tibble(rownames = "sample_id") |>
    rename(predicted = lambda.min) |>
    add_column(measured = y[-train])
  
  rsquare <- summary(lm(predicted ~ measured, data = prediction_on_test))$adj.r.squared
  
  coefs_sf <- coef(cvfit, s = "lambda.min") |>
    as.matrix() |>
    as_tibble(rownames = "transcript_id")
  
  nb_coefs <- sum(coefs_sf$s1 != 0)
  
  list(rsquare = rsquare, nb_coefs = nb_coefs,
       prediction_on_test = prediction_on_test,
       coefs_sf = coefs_sf)
}


regression_logistic_lasso <- function(my_ev, mat_sf_expression, quants){
  
  # get y data
  y <- quants[quants$event_id == my_ev, c("sample_id", "inclusion")] |>
    column_to_rownames("sample_id") |>
    as.matrix()
  
  x <- mat_sf_expression[rownames(y),]
  
  
  n <- nrow(x)
  train <- sample(n, round(.7*n))
  
  
  
  cvfit <- cv.glmnet(x[train,], y[train], family = "binomial", nfolds = 20)
  
  # Estimate on test data
  prediction_on_test <- predict(cvfit, newx = x[-train,], s = "lambda.min", type = "class") |>
    as.data.frame() |>
    as_tibble(rownames = "sample_id") |>
    rename(predicted = lambda.min) |>
    add_column(measured = y[-train])
  
  TPR <- sum(prediction_on_test$predicted == prediction_on_test$measured)/nrow(prediction_on_test)
  
  
  
  coefs_sf <- coef(cvfit, s = "lambda.min") |>
    as.matrix() |>
    as_tibble(rownames = "transcript_id")
  
  
  
  
  nb_coefs <- sum(coefs_sf$s1 != 0)
  
  list(TPR = TPR, nb_coefs = nb_coefs,
       prediction_on_test = prediction_on_test,
       coefs_sf = coefs_sf)
}


# similar to rowMax, return the value which is highest in absolute value
rowAbsMax <- function(x){
  stopifnot(all(vapply(x, is.numeric,
                       FUN.VALUE = logical(1L))))
  
  apply(x, 1, \(.row) .row[which.max(abs(.row))])
}


regression_multinomial_lasso <- function(my_ev, mat_sf_expression, quants){
  
  # get y data
  y <- quants[quants$event_id == my_ev, c("sample_id", "inclusion")] |>
    column_to_rownames("sample_id") |>
    as.matrix()
  
  x <- mat_sf_expression[rownames(y),]
  
  
  n <- nrow(x)
  train <- sample(n, round(.7*n))
  
  
  
  cvfit <- cv.glmnet(x[train,], y[train],family = "multinomial", nfolds = 20)
  
  # Estimate on test data
  prediction_on_test <- tibble(
    sample_id = rownames(x[-train,]),
    predicted = predict(cvfit, newx = x[-train,], s = "lambda.min", type = "class"),
    measured = y[-train]
  )
  
  TPR <- sum(prediction_on_test$predicted == prediction_on_test$measured)/nrow(prediction_on_test)
  
  
  
  coefs_sf <- coef(cvfit, s = "lambda.min") |>
    imap(\(.x, .i) as.matrix(.x) |>
           as_tibble(rownames = "transcript_id") |>
           rename({{.i}} := `1`))
  
  coefs_sf <- bind_cols(transcript_id = coefs_sf[[1]][["transcript_id"]],
                        map(coefs_sf, ~select(.x, -transcript_id)))
  
  coefs_sf <- bind_cols(coefs_sf,
                        s1 = rowAbsMax(select(coefs_sf, -transcript_id)))
  
  nb_coefs <- sum(coefs_sf$s1 != 0)
  
  list(TPR = TPR, nb_coefs = nb_coefs,
       prediction_on_test = prediction_on_test,
       coefs_sf = coefs_sf)
}



###

# Play with examples ----
# SE_136, SE_303 are good. 1060 particularly
(my_ev <- sample(events_to_keep_all, 1))
# my_ev <- "SE_1046"

events_coordinates |> filter(event_id== my_ev) |> pull(gene_id) |> i2s(gids)

quantifs |>
  filter(event_id == my_ev) |>
  ggplot() +
  theme_classic() +
  geom_col(aes(x = neuron_id, y = PSI, fill = neuron_id, group = sample_id),
           position = "dodge", show.legend = FALSE, color = 'white') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
  hues::scale_fill_iwanthue()

# ternary
quantifs_filtered |>
  filter(event_id == my_ev) |>
  ggplot() +
  theme_classic() +
  geom_col(aes(x = neuron_id, y = inclusion, fill = neuron_id, group = sample_id),
           position = "dodge", show.legend = FALSE, color = 'white') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
  hues::scale_fill_iwanthue()

# binary
quantifs_filtered_binary |>
  filter(event_id == my_ev) |>
  ggplot() +
  theme_classic() +
  geom_col(aes(x = neuron_id, y = inclusion, fill = neuron_id, group = sample_id),
           position = "dodge", show.legend = FALSE, color = 'white') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
  hues::scale_fill_iwanthue()



# Sparse LASSO

# Randomize!
#all values
# x[] <- sample(x)
#by sample
# x <- x[sample(nrow(x)),]


reg_lasso <- regression_lasso(my_ev, mat_sf_expression, quantifs)
reg_adalasso <- regression_adalasso(my_ev, mat_sf_expression, quantifs)
reg_loglasso <- regression_logistic_lasso(my_ev, mat_sf_expression, quantifs_filtered_binary)
reg_multilasso <- regression_multinomial_lasso(my_ev, mat_sf_expression, quantifs_filtered_ternary)

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

reg_loglasso$TPR
reg_multilasso$TPR

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


# Check against putative SFs



xx <- reg_lasso$coefs_sf |>
  filter(transcript_id != "(Intercept)") |>
  mutate(gene_id = convert_sf_tx2g(transcript_id),
         putative_sf = gene_id %in% putative_splice_factors)

ggplot(xx) +
  theme_classic() +
  # geom_boxplot(aes(x = putative_sf, y = s1))
  geom_jitter(aes(x = putative_sf, y = s1), alpha = .2)

xx |>
  filter(s1 != 0) |>
  ggplot() +
  theme_classic() +
  geom_jitter(aes(x = putative_sf, y = s1))

eulerr::euler(list(lasso = unique(xx$gene_id[xx$s1 != 0]),
                   literature = putative_splice_factors)) |>
  plot(quantities = TRUE)



## Continue
reg_loglasso$coefs_sf |>
  filter(s1 != 0) |>
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
reg_multilasso$coefs_sf |>
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
     logistic = reg_loglasso$coefs_sf |> filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     multinomial = reg_multilasso$coefs_sf |> filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     literature = putative_splice_factors) |>
  eulerr::euler() |>
  plot(quantities = TRUE)



list(lasso    = reg_lasso$coefs_sf |>  filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     adalasso = reg_adalasso$coefs_sf |> filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     logistic = reg_loglasso$coefs_sf |> filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     multinomial = reg_multilasso$coefs_sf |> filter(s1 != 0) |> pull(transcript_id) |> setdiff("(Intercept)"),
     literature = putative_splice_factors) |>
  UpSetR::fromList() |>
  UpSetR::upset()



# First pass ----
# make a first round of regressions


first_pass_regression <- quantifs_filtered |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  arrange(desc(nb_samples), desc(nb_neurons)) |>
  # slice_head(n = 2) |>
  mutate(fit_lasso = map(event_id,
                         \(ev) possibly(regression_lasso,
                                        otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                         .progress = TRUE),
         fit_adalasso = map(event_id,
                            \(ev) possibly(regression_adalasso,
                                           otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                            .progress = TRUE),
         fit_loglasso = map(event_id,
                            \(ev) possibly(regression_logistic_lasso,
                                           otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered_binary),
                            .progress = TRUE),
         fit_multlasso = map(event_id,
                             \(ev) possibly(regression_multinomial_lasso,
                                            otherwise = list(NA))(ev, mat_sf_expression, quantifs_filtered),
                             .progress = TRUE))
# qs::qsave(first_pass_regression, "data/intermediates/230130_fits_for_quantifs_4methods_allgenes_cache.qs")


first_pass_regression <- qs::qread("data/intermediates/230125_fits_for_quantifs_4methods_cache.qs")


# Compare results on first pass ----

patchwork::wrap_plots(
  ggplot(first_pass_regression |>
           mutate(rsquare = map_dbl(fit_lasso, pluck("rsquare"))),
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(nb_coefs = map_int(fit_lasso, pluck("nb_coefs"))),
         aes(x = nb_samples, y = nb_coefs)) +
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
           mutate(nb_coefs = map_int(fit_adalasso, pluck("nb_coefs"))),
         aes(x = nb_samples, y = nb_coefs)) +
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
           mutate(nb_coefs = map_int(fit_loglasso, ~ pluck(., "nb_coefs", .default = NA_integer_))),
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression |>
           mutate(nb_coefs = map_int(fit_multlasso, ~ pluck(., "nb_coefs", .default = NA_integer_))),
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ncol = 2
)




first_pass_regression |>
  mutate(lasso_rsquare = map_dbl(fit_lasso, ~ pluck(., "rsquare", .default = NA_integer_)),
         adalasso_rsquare = map_dbl(fit_adalasso, ~ pluck(., "rsquare", .default = NA_integer_))) |>
  select(event_id, ends_with("_rsquare")) |>
  pivot_longer(-event_id, names_to = "method", values_to = "rsquare", names_pattern = "(a?d?a?lasso)_rsquare") |>
  ggplot(aes(x = method, y = rsquare)) +
  theme_classic() +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm()

xx <- first_pass_regression |>
  mutate(lasso_coefs_sf = map(fit_lasso, ~ pluck(., "coefs_sf", .default = tibble())),
         adalasso_coefs_sf = map(fit_adalasso, ~ pluck(., "coefs_sf", .default = tibble())),
         loglasso_coefs_sf = map(fit_loglasso, ~ pluck(., "coefs_sf", .default = tibble())),
         multlasso_coefs_sf = map(fit_multlasso, ~ pluck(., "coefs_sf", .default = tibble()))) |>
  select(event_id, ends_with("_coefs_sf")) |>
  pivot_longer(-event_id, names_to = "method", values_to = "coefs", names_pattern = "([admultog]{0,4}lasso)_coefs_sf") |>
  unnest(coefs) |>
  filter(s1 != 0,
         transcript_id != "(Intercept)")
  
xx |>
  ggplot(aes(x = method, y = abs(s1))) +
  theme_classic() +
  geom_boxplot() +
  scale_y_log10()

ggplot(xx, aes(x = transcript_id, y = s1)) +
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








