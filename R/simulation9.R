# Simulation to check whether the permutation test gives results that make sense

# Version 9, expand on the previous versions (`simulation2.R`) for a more realistic simulated dataset

## --- Approach ---
# Keep measured TPMs, (no randomization)
# For each event, select 3+3 transcripts non-zero, all other coefficients 0
# For these 6 “true non-0” transcripts, coefficient from Unif(-5, 5)
# Compute simulation:
#  number of included reads and number of excluded reads simulated from 3 SFs
#  Nincl = NB(mu, size) where mu = m0 + coef*log(TPM)
#  Nexcl = NB(mu, size) where mu = m0 - coef*log(TPM)
# where m0 same for incl and excl, represents TF-controlled gene expression level
# Also check size and m0 relationship
#  then PSI = Nincl/(Nincl+Nexcl)
# 
# Create 100 simulated datasets
# 










# ----Simulations ----

library(tidyverse)

# read real data to generate simulations from
# we use the pre-filtered events from `prefilter_quantifs.R`

events_to_keep <- read_lines("data/intermediates/230512_simulation_v10/events_to_keep.txt")
neurons_to_keep <- read_lines("data/intermediates/230512_simulation_v10/neurons_to_keep.txt")

quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) |>
  filter(! is.na(PSI),
         event_id %in% events_to_keep)

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") |>
  filter(gene_id %in% putative_splice_factors) |>
  filter(sample_id %in% unique(quantifs$sample_id)) # remove RICr133, Ref, ...




# TPM heatmap
mat_sf_expression <- sf_expression |>
  pivot_wider(id_cols = transcript_id,
              names_from = sample_id,
              values_from = TPM) |>
  column_to_rownames("transcript_id") |>
  as.matrix()

pheatmap::pheatmap(log1p(mat_sf_expression),scale = "column",
                   annotation_col = sf_expression |> select(sample_id, neuron_id) |> distinct() |> column_to_rownames("sample_id"))




# Parameters of real data ----

counts_measured <- quantifs |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads)) |>
  select(event_id, sample_id, Nincl, Nexcl) |>
  pivot_longer(c(Nincl, Nexcl), names_to = "contribution", values_to = "count") |>
  group_by(event_id, contribution) |>
  nest()

# fit NB to measured counts
all_fits_res <- map(counts_measured$data,
           ~ possibly(\(dat) fitdistrplus::fitdist(data = dat[["count"]],
                                                   distr = "nbinom"))(.x)) |>
  set_names(paste0(counts_measured$event_id,"_",counts_measured$contribution)) |>
  rlist::list.clean() |>
  imap(~ tibble(name = .y, mu = .x$estimate[["mu"]], size = .x$estimate[["size"]])) |>
  list_rbind() |>
  extract(name,
          into = c("event_id", "contribution"),
          regex = "^(SE_[[:digit:]]{1,4})_(Nincl|Nexcl)$")


# With separate Nincl and Nexcl

all_fits_res |>
  pivot_longer(c(mu, size), names_to = "parameter",values_to = "value") |>
  ggplot() + theme_classic() +
  geom_density(aes(x = value, fill = contribution), alpha = .2) +
  scale_x_log10() +
  facet_wrap(~parameter, scales = "free_x") +
  geom_vline(aes(xintercept = median),
             data = tibble(parameter = c("mu", "size"),
                           median = c(median(all_fits_res$mu), median(all_fits_res$size))))


# fit log-normal to mu
fit_mu_incl <- fitdistrplus::fitdist(data = all_fits_res |> filter(contribution == "Nincl") |> pull(mu),
                                     distr = "lnorm")
plot(fit_mu_incl, breaks = 50)



fit_mu_excl <- fitdistrplus::fitdist(data = all_fits_res |> filter(contribution == "Nexcl") |> pull(mu),
                                     distr = "lnorm")
plot(fit_mu_excl, breaks = 50)


# Parameters
#   mu     meanlog       sdlog
# incl    3.533800    1.004089
# excl    3.512100    1.047392 

# check mu0
all_fits_res |>
  select(-size) |>
  pivot_wider(id_cols = event_id,
              names_from = "contribution",
              values_from = "mu") |>
  filter(Nincl < 1000 & Nexcl < 1000) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = Nincl, y = Nexcl))


# compute mu0
all_fits_mu0 <- all_fits_res |>
  select(-size) |>
  pivot_wider(id_cols = event_id,
              names_from = "contribution",
              values_from = "mu") |>
  mutate(mu0 = (Nincl + Nexcl)/2,
         S = (Nincl - Nexcl)/(Nincl + Nexcl)) |>
  right_join(all_fits_res, by = "event_id")

hist(log(all_fits_mu0$mu0), breaks = 50)
hist(all_fits_mu0$S, breaks = 50)
plot(log(all_fits_mu0$mu0), all_fits_mu0$S)

fit_mu0 <- fitdistrplus::fitdist(data = all_fits_mu0$mu0,
                                 distr = "lnorm")
plot(fit_mu0, breaks = 50)

fit_s <- fitdistrplus::fitdist(data = all_fits_mu0$S,
                               distr = "unif")
plot(fit_s, breaks = 50)



# Size: relationship vs mu0 and S
all_fits_expanded <- all_fits_res |>
  left_join(all_fits_mu0 |> select(event_id, mu0, S),
            by = "event_id",
            relationship = "many-to-many") 

patchwork::wrap_plots(
  all_fits_mu0 |>
  filter(mu < 1000) |>
  ggplot(aes(x = mu0, y = size, color = contribution)) +
  theme_classic() +
  geom_point(show.legend = FALSE) +
  scale_x_log10() +
  geom_smooth(method = lm, formula = y ~ x, show.legend = FALSE),
  all_fits_mu0 |>
    filter(mu < 1000) |>
    ggplot(aes(x = S, y = size, color = contribution)) +
    theme_classic() +
    geom_point() +
    geom_smooth(method = lm, formula = y ~ x)
)


# fits size to predict it from S
mod_size_mu_incl <- lm(size ~ S, data = all_fits_mu0 |> filter(contribution == "Nincl"))
mod_size_mu_excl <- lm(size ~ S, data = all_fits_mu0 |> filter(contribution == "Nexcl"))

plot(mod_size_mu_incl)
plot(mod_size_mu_excl)


# simulate from the fits to check (function defined below)
all_fits_res |> select(-event_id) |> add_column(type = "real") |>
  bind_rows(tibble(type = "simulated",
                   contribution = "Nincl",
                   size = simul_size(all_fits_mu0$S[all_fits_mu0$contribution == "Nincl"], mod_size_mu_incl))) |>
  bind_rows(tibble(type = "simulated",
                   contribution = "Nexcl",
                   size = simul_size(all_fits_mu0$S[all_fits_mu0$contribution == "Nexcl"], mod_size_mu_excl))) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~contribution) +
  # geom_freqpoly(aes(x = size, color = type))
  geom_density(aes(x = size, fill = type), alpha = .2, bw = .1)



real_data_fit <- list(mu0 = fit_mu0$estimate,
                      mod_size_mu_incl = mod_size_mu_incl,
                      mod_size_mu_excl = mod_size_mu_excl)
# qs::qsave(real_data_fit,
#           "data/intermediates/230512_simulation_v10/real_data_fit.qs")







# Simulate single ----

library(tidyverse)

# read real data to generate simulations from
# we use the pre-filtered data from `prefilter_quantifs.R`

events_to_keep <- read_lines("data/intermediates/230512_simulation_v10/events_to_keep.txt")
neurons_to_keep <- read_lines("data/intermediates/230512_simulation_v10/neurons_to_keep.txt")

quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) |>
  filter(! is.na(PSI),
         event_id %in% events_to_keep)

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") |>
  filter(gene_id %in% putative_splice_factors) |>
  filter(sample_id %in% unique(quantifs$sample_id)) # remove RICr133, Ref, ...

# prepare simulation data
real_data_fit <- qs::qread("data/intermediates/230512_simulation_v10/real_data_fit.qs")


#~ functions ----

# copied from `simulate()`: given mu, predict (with realistic noise) a corresponding size.
simul_size <- function(mu, fit){
  vars <- deviance(fit)/df.residual(fit)
  sim <- coef(fit)["(Intercept)"] + coef(fit)["S"] * mu + rnorm(length(mu), sd = sqrt(vars))
  # censor negative values
  if(any(sim <= 0)) sim[sim <= 0] <- simul_size(mu[sim <= 0], fit)
  sim
}

rcoefs <- function(n, nb_nonzero = 3L, range_unif = 5, noise = 0.02){
  stopifnot(range_unif > 0)
  unif_coefs <- runif(n, min = -range_unif, max = range_unif)
  noise_coefs <- rnorm(n, mean = 0, sd = noise)
  non_zeros <- rep(c(1, 0), times =  c(nb_nonzero, n - nb_nonzero)) |> sample()
  
  non_zeros * unif_coefs
}

# rescale a centered normal into a uniform
# see https://math.stackexchange.com/questions/1063865/transforming-a-normal-distribution-to-a-uniform-one
# and https://math.stackexchange.com/questions/2343952/how-to-transform-gaussiannormal-distribution-to-uniform-distribution
rescale_distr <- function(x){
  2*pnorm(x/sd(x)) - 1
}


#~ descriptors ----
nb_tx  <- sf_expression |>
  pull(transcript_id) |>
  unique() |>
  length()
nb_events <- quantifs |>
  pull(event_id) |>
  unique() |>
  length()
nb_samples <- quantifs |>
  pull(sample_id) |>
  unique() |>
  length()
nb_datapoints <- quantifs |> nrow()

sim_sf <- sf_expression |>
  select(transcript_id, sample_id, neuron_id, TPM)

true_coefs <- expand_grid(event_id = unique(quantifs$event_id),
                          transcript_id = unique(sim_sf$transcript_id)) |>
  group_by(event_id) |>
  nest() |>
  mutate(data = map(data, ~ add_column(.x, true_coef = rcoefs(nb_tx, nb_nonzero = 40, range_unif = 50, noise = 0)))) |>
  unnest(data) |> ungroup()

# get list of coefficients for each sample, event, transcript; generate mu0
randomized_samples <- quantifs |>
  select(event_id, sample_id) |>
  mutate(mu0 = rlnorm(n = nrow(quantifs),
                      meanlog = real_data_fit$mu0[["meanlog"]],
                      sdlog = real_data_fit$mu0[["sdlog"]])) |>
  full_join(sim_sf, by = "sample_id", relationship = "many-to-many") |>
  left_join(true_coefs, by = c("event_id", "transcript_id"), relationship = "many-to-many")


# check distrib S as simulated
sim_quantifs_raw <- randomized_samples |>
  group_by(event_id, sample_id, mu0) |>
  summarize(Sraw = sum((true_coef*log1p(TPM))),
            .groups = 'drop')

hist(sim_quantifs_raw$Sraw, breaks = 50)
hist(rescale_distr(sim_quantifs_raw$Sraw), breaks = 50)




# compute
sim_quantifs <- sim_quantifs_raw |>
  mutate(S = rescale_distr(Sraw),
         mu_incl = mu0 * (1 + S),
         mu_excl = mu0 * (1 - S),
         size_incl = simul_size(S, real_data_fit$mod_size_mu_incl),
         size_excl = simul_size(S, real_data_fit$mod_size_mu_excl),
         N_incl = rnbinom(n = nb_datapoints,
                         size = size_incl,
                         mu = mu_incl),
         N_excl = rnbinom(n = nb_datapoints,
                         size = size_excl,
                         mu = mu_excl)) |>
  mutate(PSI = N_incl / (N_incl + N_excl),
         nb_reads = N_incl + N_excl,
         neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{2,4}$")[,2]) |>
  select(event_id, sample_id, nb_reads, PSI, neuron_id)


# Compare histograms with real data
tibble(type = c(rep("measured", nb_datapoints), rep("simulated", nb_datapoints)),
       PSI = c(quantifs$PSI, sim_quantifs$PSI)) |>
  ggplot() + theme_classic() +
  geom_density(aes(x = PSI, fill = type), alpha = .3, bw = .01)

tibble(type = c(rep("measured", nb_datapoints), rep("simulated", nb_datapoints)),
       PSI = c(quantifs$PSI, sim_quantifs$PSI)) |>
  ggplot() + theme_classic() +
  geom_freqpoly(aes(x = PSI, color = type), bins = 100)


tibble(type = c(rep("measured", nb_datapoints), rep("simulated", nb_datapoints)),
       `Total number of reads` = c(quantifs$nb_reads, sim_quantifs$nb_reads)) |>
  ggplot() + theme_classic() +
  geom_density(aes(x = `Total number of reads`, fill = type), alpha = .3, bw = .01) +
  scale_x_log10()

tibble(type = c(rep("measured", nb_datapoints), rep("simul", nb_datapoints)),
       `Total number of reads` = c(quantifs$nb_reads, sim_quantifs$nb_reads)) |>
  ggplot() + theme_classic() +
  geom_freqpoly(aes(x = `Total number of reads`, color = type), bins = 100) +
  scale_x_log10()



# Filter ----

#~ Filter neurons ----

# Filter neurons with too few samples. Also remove Ref
keep_neurons_real <- sf_expression |>
  select(sample_id, neuron_id) |>
  distinct() |>
  dplyr::count(neuron_id) |>
  filter(n > 2) |>
  pull(neuron_id) |>
  setdiff("Ref")

keep_neurons_sim <- sim_sf |>
  select(sample_id, neuron_id) |>
  distinct() |>
  dplyr::count(neuron_id) |>
  filter(n > 2) |>
  pull(neuron_id)




sf_expression <- sf_expression |>
  filter(neuron_id %in% keep_neurons_real)
sim_sf <- sim_sf |>
  filter(neuron_id %in% keep_neurons_sim)


quantifs <- quantifs |>
  filter(neuron_id %in% keep_neurons_real)
sim_quantifs <- sim_quantifs |>
  filter(neuron_id %in% keep_neurons_sim)



#~ Filter events ----

# Keep only if enough reads supporting measure
quantifs |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads),color='grey', bins = 200) +
  geom_vline(aes(xintercept = 20), linetype = 'dashed') +
  scale_x_log10()

quantifs_filtered_n_reads_real <- quantifs |>
  filter(nb_reads > 20)


sim_quantifs |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads),color='grey', bins = 200) +
  geom_vline(aes(xintercept = 20), linetype = 'dashed') +
  scale_x_log10()

quantifs_filtered_n_reads_sim <- sim_quantifs |>
  filter(nb_reads > 20)




# filter based on nb of samples that event was measured in
quantifs_filtered_n_reads_real |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  ggplot() + theme_classic() +
  geom_point(aes(x = nb_neurons, y = nb_samples), alpha = .2) +
  xlab("Number of neurons") + ylab("Number of samples") +
  geom_hline(aes(yintercept = 100.5), color = 'darkred') +
  geom_vline(aes(xintercept = 32.5), color = 'darkred')



events_to_keep_n_samples_real <- quantifs_filtered_n_reads_real |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  filter(nb_samples > 100,
         nb_neurons > 32) |>
  pull(event_id)


quantifs_filtered_nsamples_real <- quantifs_filtered_n_reads_real |>
  filter(event_id %in% events_to_keep_n_samples_real)








quantifs_filtered_n_reads_sim |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  ggplot() + theme_classic() +
  geom_point(aes(x = nb_neurons, y = nb_samples), alpha = .2) +
  xlab("Number of neurons") + ylab("Number of samples") +
  geom_hline(aes(yintercept = 100.5), color = 'darkred') +
  geom_vline(aes(xintercept = 32.5), color = 'darkred')



events_to_keep_n_samples_sim <- quantifs_filtered_n_reads_sim |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  filter(nb_samples > 100,
         nb_neurons > 32) |>
  pull(event_id)


quantifs_filtered_nsamples_sim <- quantifs_filtered_n_reads_sim |>
  filter(event_id %in% events_to_keep_n_samples_sim)



# Filter events to remove those that are not DS between neuron types
quantifs_filtered_nsamples_real |>
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
  geom_hline(aes(yintercept = 0.05), color = 'grey') +
  ggtitle("Real data")

events_to_keep_variability_real <- quantifs_filtered_nsamples_real |>
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


quantifs_filtered_real <- quantifs_filtered_nsamples_real |>
  filter(event_id %in% events_to_keep_variability_real) |>
  filter(! is.na(PSI))

quantifs_filtered_real |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads), color='grey', bins = 100) +
  scale_x_log10()







quantifs_filtered_nsamples_sim |>
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
  geom_hline(aes(yintercept = 0.05), color = 'grey') +
  ggtitle("Simulated data")

events_to_keep_variability_sim <- quantifs_filtered_nsamples_sim |>
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


quantifs_filtered_sim <- quantifs_filtered_nsamples_sim |>
  filter(event_id %in% events_to_keep_variability_sim) |>
  filter(! is.na(PSI))



quantifs_filtered_sim |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads),color='grey', bins = 100) +
  scale_x_log10()






#~ Check match after filtering ----
nb_datapoints_real <- nrow(quantifs_filtered_real)
nb_datapoints_sim <- nrow(quantifs_filtered_sim)

tibble(type = c(rep("measured", nb_datapoints_real), rep("simul", nb_datapoints_sim)),
       PSI = c(quantifs_filtered_real$PSI, quantifs_filtered_sim$PSI)) |>
  ggplot() + theme_classic() +
  geom_density(aes(x = PSI, fill = type), alpha = .3, bw = .01)

tibble(type = c(rep("measured", nb_datapoints_real), rep("simul", nb_datapoints_sim)),
       PSI = c(quantifs_filtered_real$PSI, quantifs_filtered_sim$PSI)) |>
  ggplot() + theme_classic() +
  geom_freqpoly(aes(x = PSI, color = type), bins = 100)


tibble(type = c(rep("measured", nb_datapoints_real), rep("simul", nb_datapoints_sim)),
       `Total number of reads` = c(quantifs_filtered_real$nb_reads, quantifs_filtered_sim$nb_reads)) |>
  ggplot() + theme_classic() +
  geom_density(aes(x = `Total number of reads`, fill = type), alpha = .3, bw = .01) +
  scale_x_log10()

tibble(type = c(rep("measured", nb_datapoints_real), rep("simul", nb_datapoints_sim)),
       `Total number of reads` = c(quantifs_filtered_real$nb_reads, quantifs_filtered_sim$nb_reads)) |>
  ggplot() + theme_classic() +
  geom_freqpoly(aes(x = `Total number of reads`, color = type), bins = 100) +
  scale_x_log10()







qs::qsave(quantifs_filtered_sim, "data/intermediates/230517_simulation_v12/quantifs_filtered.qs")
qs::qsave(sim_sf, "data/intermediates/230517_simulation_v12/sim_sf.qs")
qs::qsave(true_coefs, "data/intermediates/230517_simulation_v12/true_coefs.qs")




# Regression on single simul ----

library(tidyverse)
library(glmnet)

source("R/regression_functions.R")


# Read data ----
sim_quantifs <- qs::qread("data/intermediates/230517_simulation_v12/quantifs_filtered.qs")
sim_sf <- qs::qread("data/intermediates/230517_simulation_v12/sim_sf.qs")
sim_true_coefs <- qs::qread("data/intermediates/230517_simulation_v12/true_coefs.qs")



# Predict DeltaPSI instead of PSI
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9 for rationale

logit <- function(x){
  stopifnot(all((x >= 0 & x <= 1) | is.nan(x)))
  x[x == 0] <- .01
  x[x == 1] <- .99
  log(x/(1-x))
}

sim_quantifs <- sim_quantifs |>
                      group_by(event_id) |>
                      mutate(dPSI_nat = PSI - mean(PSI, na.rm = TRUE),
                             dPSI_logit = logit(PSI) - logit(mean(PSI, na.rm = TRUE))) |>
  ungroup()


events_to_keep <- unique(sim_quantifs$event_id)





# Make SF expression as a matrix for use in regression

mat_sf_expression <-  sim_sf |>
  select(transcript_id, sample_id, TPM) |>
  mutate(TPM = log(TPM + 1)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "TPM") |>
  column_to_rownames("sample_id") |>
  as.matrix() |>
  scale()
mat_sf_expression <- mat_sf_expression[,! apply(mat_sf_expression, 2, \(col) any(is.na(col)))]



# Play with examples ----
# SE_136, SE_303 are good. 1060 particularly
(my_ev <- sample(events_to_keep, 1))
# my_ev <- "SE_1046"


sim_quantifs |>
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
y <- sim_quantifs[sim_quantifs$event_id == my_ev, c("sample_id", "dPSI_nat")] |>
  column_to_rownames("sample_id") |>
  filter(!is.na(dPSI_nat)) |>
  as.matrix()

# get x data
x <- mat_sf_expression[rownames(y),]

n <- nrow(x)
train <- sample(n, round(.7*n))

fit <- cv.glmnet(x = x[train,], y = y[train,], type.measure = "mse", nfolds = 20,
                 alpha = 1, intercept = FALSE, standardize = TRUE)
plot(fit)


# check on test set
prediction_on_test <- predict(fit, newx = x[-train,], s = "lambda.min") |>
  as.data.frame() |>
  as_tibble(rownames = "sample_id") |>
  rename(predicted = lambda.min) |>
  add_column(measured = y[-train])

prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")

coefs_sf <- coef(fit, s = "lambda.min") |>
  as.matrix() |>
  as_tibble(rownames = "transcript_id")

coefs_sf |>
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

true_coefs <- sim_true_coefs |>
  filter(event_id == my_ev) |>
  group_by(event_id, transcript_id) |>
  summarize(true_coef = sum(true_coef),
            .groups = 'drop')


full_join(coefs_sf |> filter(transcript_id != "(Intercept)") |> select(transcript_id, computed = s1),
          true_coefs |> select(transcript_id, ground_truth = true_coef),
          by = "transcript_id") |>
  ggplot() + theme_classic() +
  geom_vline(aes(xintercept = 0), color = 'grey') +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = ground_truth, y = computed))


# check influence visually (i.e. check that simulation not too noisy)
true_coefs |> arrange(desc(abs(true_coef))) |> head(3)



left_join(
  left_join(
    sim_sf |>
      filter(transcript_id == "STRG.2883.2") |>
      select(sample_id, TPM_D = TPM),
    sim_quantifs |>
      ungroup() |>
      filter(event_id == my_ev) |>
      select(sample_id, PSI),
    by = "sample_id"
  ),
  left_join(
    sim_sf |>
      filter(transcript_id == "STRG.1415.4") |>
      select(sample_id, TPM_C = TPM),
    sim_quantifs |>
      ungroup() |>
      filter(event_id == my_ev) |>
      select(sample_id, PSI),
    by = "sample_id"
  ),
  by = c("sample_id", "PSI")
) |>
  plotly::plot_ly(x = ~TPM_D, y = ~TPM_C, z = ~PSI)



# try LASSO approaches
(my_ev <- sample(events_to_keep, 1))
reg_lasso1 <- regression_wrapper(my_ev = my_ev,
                                 regression_method = "lasso",
                                 column = "PSI",
                                 shuffle = FALSE,
                                 mat_sf_expression = mat_sf_expression,
                                 quants = sim_quantifs,
                                 intercept = TRUE)

reg_lasso1$prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")

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
true_coefs1 <- sim_true_coefs |>
  filter(event_id == my_ev) |>
  group_by(event_id, transcript_id) |>
  summarize(true_coef = sum(true_coef),
            .groups = 'drop')
# mutate(true_coef = if_else(abs(true_coef) > 0.8,
#                            true_coef, 0))

true_coefs1[true_coefs1$true_coef != 0,]


















# Simulate 100 datasets ----


library(tidyverse)

events_to_keep <- read_lines("data/intermediates/230512_simulation_v10/events_to_keep.txt")
neurons_to_keep <- read_lines("data/intermediates/230512_simulation_v10/neurons_to_keep.txt")

quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) |>
  filter(! is.na(PSI),
         event_id %in% events_to_keep)

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") |>
  filter(gene_id %in% putative_splice_factors) |>
  filter(sample_id %in% unique(quantifs$sample_id)) # remove RICr133, Ref, ...

real_data_fit <- qs::qread("data/intermediates/230512_simulation_v10/real_data_fit.qs")







# copied from `simulate()`: given mu, predict (with realistic noise) a corresponding size.
simul_size <- function(mu, fit){
  vars <- deviance(fit)/df.residual(fit)
  sim <- coef(fit)["(Intercept)"] + coef(fit)["S"] * mu + rnorm(length(mu), sd = sqrt(vars))
  # censor negative values
  if(any(sim <= 0)) sim[sim <= 0] <- simul_size(mu[sim <= 0], fit)
  sim
}

rcoefs <- function(n, nb_nonzero = 3L, range_unif = 5, noise = 0.02){
  stopifnot(range_unif > 0)
  unif_coefs <- runif(n, min = -range_unif, max = range_unif)
  noise_coefs <- rnorm(n, mean = 0, sd = noise)
  non_zeros <- rep(c(1, 0), times =  c(nb_nonzero, n - nb_nonzero)) |> sample()
  
  non_zeros * unif_coefs
}

# rescale a centered normal into a uniform
# see https://math.stackexchange.com/questions/1063865/transforming-a-normal-distribution-to-a-uniform-one
# and https://math.stackexchange.com/questions/2343952/how-to-transform-gaussiannormal-distribution-to-uniform-distribution
rescale_distr <- function(x){
  2*pnorm(x/sd(x)) - 1
}



simulate_single <- function(real_data_fit, sf_expression, quantifs){
  
  #~ descriptors ----
  nb_tx  <- sf_expression |>
    pull(transcript_id) |>
    unique() |>
    length()
  
  nb_datapoints <- quantifs |> nrow()
  
  sim_sf <- sf_expression |>
    select(transcript_id, sample_id, neuron_id, TPM)
  
  
  true_coefs <- expand_grid(event_id = unique(quantifs$event_id),
                            transcript_id = unique(sim_sf$transcript_id)) |>
    group_by(event_id) |>
    nest() |>
    mutate(data = map(data, ~ add_column(.x, true_coef = rcoefs(nb_tx, nb_nonzero = 40)))) |>
    unnest(data) |> ungroup()
  
  # get list of coefficients for each sample, event, transcript; generate mu0
  randomized_samples <- quantifs |>
    select(event_id, sample_id) |>
    mutate(mu0 = rlnorm(n = nrow(quantifs),
                        meanlog = real_data_fit$mu0[["meanlog"]],
                        sdlog = real_data_fit$mu0[["sdlog"]])) |>
    full_join(sim_sf, by = "sample_id", relationship = "many-to-many") |>
    left_join(true_coefs, by = c("event_id", "transcript_id"), relationship = "many-to-many")
  
  # check distrib S as simulated
  sim_quantifs_raw <- randomized_samples |>
    group_by(event_id, sample_id, mu0) |>
    summarize(Sraw = sum((true_coef*log1p(TPM))),
              .groups = 'drop')
  
  # compute
  sim_quantifs <- sim_quantifs_raw |>
    mutate(S = rescale_distr(Sraw),
           mu_incl = mu0 * (1 + S),
           mu_excl = mu0 * (1 - S),
           size_incl = simul_size(S, real_data_fit$mod_size_mu_incl),
           size_excl = simul_size(S, real_data_fit$mod_size_mu_excl),
           N_incl = rnbinom(n = nb_datapoints,
                            size = size_incl,
                            mu = mu_incl),
           N_excl = rnbinom(n = nb_datapoints,
                            size = size_excl,
                            mu = mu_excl)) |>
    mutate(PSI = N_incl / (N_incl + N_excl),
           nb_reads = N_incl + N_excl,
           neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{2,4}$")[,2]) |>
    select(event_id, sample_id, nb_reads, PSI, neuron_id)
  
  list(sim_quantifs = sim_quantifs, sim_sf = sim_sf, true_coefs = true_coefs)
}

# 20s per replicate, 30' for 100
sim_replicated <- map(1:100,
                      ~simulate_single(real_data_fit, sf_expression, quantifs),
                      .progress = TRUE)


qs::qsave(sim_replicated, "data/intermediates/230517_simulation_v11/rep_simulations.qs")






# ---- Regression ----

library(tidyverse)
library(glmnet)

source("R/regression_functions.R")


filter_single <- function(sims){
  
  stopifnot(class(sims) == "list")
  stopifnot(length(sims) == 3L)
  stopifnot(all.equal(names(sims), c("sim_quantifs", "sim_sf", "true_coefs")))
  
  names(sims) <- c("quantifs", "sf", "coefs")
  
  ## Filter neurons
  keep_neurons <- sims$sf |>
    select(sample_id, neuron_id) |>
    distinct() |>
    dplyr::count(neuron_id) |>
    filter(n > 2) |>
    pull(neuron_id)
  
  sims$sf <- sims$sf |>
    filter(neuron_id %in% keep_neurons)
  
  sims$quantifs <- sims$quantifs |>
    filter(neuron_id %in% keep_neurons)
  
  
  ## Filter events
  
  # Keep only if enough reads supporting measure
  sims$quantifs <- sims$quantifs |>
    filter(nb_reads > 20)
  
  # filter based on nb of samples that event was measured in
  events_to_keep_n_samples <- sims$quantifs |>
    filter(! is.na(PSI)) |>
    group_by(event_id) |>
    summarize(nb_samples = n(),
              nb_neurons = n_distinct(neuron_id)) |>
    filter(nb_samples > 100,
           nb_neurons > 32) |>
    pull(event_id)
  
  sims$quantifs <- sims$quantifs |>
    filter(event_id %in% events_to_keep_n_samples)
  
  # Filter events to remove those that are not DS between neuron types
  events_to_keep_variability <- sims$quantifs |>
    filter(!is.na(PSI)) |>
    group_by(neuron_id, event_id) |>
    summarize(mean_PSI  = mean(PSI),
              sd_PSI = sd(PSI),
              .groups = 'drop') |>
    group_by(event_id) |>
    summarize(mean_PSI_btw_neurs  = mean(mean_PSI),
              sd_PSI_btw_neurs = sd(mean_PSI),
              .groups = 'drop') |>
    filter(sd_PSI_btw_neurs > 0.05) |>
    pull(event_id)
  
  sims$quantifs <- sims$quantifs |>
    filter(event_id %in% events_to_keep_variability) |>
    filter(! is.na(PSI))
  
  list(sim_quantifs = sims$quantifs, sim_sf = sims$sf, true_coefs = sims$coefs)
}



#~ Read data ----
sim_replicated <- qs::qread("data/intermediates/230512_simulation_v10/rep_simulations.qs")


sim_filtered <- sim_replicated |>
  map(filter_single,
      .progress = TRUE)


sim_quantifs <- list_transpose(sim_filtered)[["sim_quantifs"]]
sim_sf <- list_transpose(sim_filtered)[["sim_sf"]]
sim_true_coefs <- list_transpose(sim_filtered)[["true_coefs"]]



# Predict DeltaPSI instead of PSI
# see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02273-7#Sec9 for rationale

logit <- function(x){
  stopifnot(all((x >= 0 & x <= 1) | is.nan(x)))
  x[x == 0] <- .01
  x[x == 1] <- .99
  log(x/(1-x))
}

sim_quantifs <- map(sim_quantifs,
                    ~ .x |>
                      group_by(event_id) |>
                      mutate(dPSI_nat = PSI - mean(PSI, na.rm = TRUE),
                             dPSI_logit = logit(PSI) - logit(mean(PSI, na.rm = TRUE))))


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
                  mat_sf_expression[,! apply(mat_sf_expression, 2, \(col) any(is.na(col)))]},
                  .progress = TRUE)



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
y <- quantifs_filtered[quantifs_filtered$event_id == my_ev, c("sample_id", "dPSI_nat")] |>
  column_to_rownames("sample_id") |>
  filter(!is.na(dPSI_nat)) |>
  as.matrix()

# get x data
x <- mat_sf_expression[rownames(y),]

n <- nrow(x)
train <- sample(n, round(.7*n))

fit <- cv.glmnet(x = x[train,], y = y[train,], type.measure = "mse", nfolds = 20, alpha = 1, intercept = FALSE, standardize = TRUE)
plot(fit)


# check on test set
prediction_on_test <- predict(fit, newx = x[-train,], s = "lambda.min") |>
  as.data.frame() |>
  as_tibble(rownames = "sample_id") |>
  rename(predicted = lambda.min) |>
  add_column(measured = y[-train])

prediction_on_test |>
  ggplot(aes(x = measured, y = predicted)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm")

coefs_sf <- coef(fit, s = "lambda.min") |>
  as.matrix() |>
  as_tibble(rownames = "transcript_id")

coefs_sf |>
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

true_coefs <- sim_true_coefs[[my_rep]] |>
  filter(event_id == my_ev) |>
  group_by(event_id, transcript_id) |>
  summarize(true_coef = sum(true_coef),
            .groups = 'drop')


full_join(coefs_sf |> filter(transcript_id != "(Intercept)") |> select(transcript_id, computed = s1),
          true_coefs |> select(transcript_id, ground_truth = true_coef),
          by = "transcript_id") |>
  ggplot() + theme_classic() +
  geom_vline(aes(xintercept = 0), color = 'grey') +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_point(aes(x = ground_truth, y = computed))


# check influence visually (i.e. check that simulation not too noisy)
true_coefs |> arrange(desc(abs(true_coef)))


left_join(
  sim_sf[[my_rep]] |>
    filter(transcript_id == "K02H8.1d.1") |>
    select(sample_id, TPM),
  sim_quantifs[[my_rep]] |>
    ungroup() |>
    filter(event_id == my_ev) |>
    select(sample_id, PSI),
  by = "sample_id"
) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = log(TPM), y = PSI))

left_join(
  sim_sf[[my_rep]] |>
    filter(transcript_id == "STRG.6313.3") |>
    select(sample_id, TPM),
  sim_quantifs[[my_rep]] |>
    ungroup() |>
    filter(event_id == my_ev) |>
    select(sample_id, PSI),
  by = "sample_id"
) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = log(TPM), y = PSI))


left_join(
  left_join(
    sim_sf[[my_rep]] |>
      filter(transcript_id == "K02H8.1d.1") |>
      select(sample_id, TPM_D = TPM),
    sim_quantifs[[my_rep]] |>
      ungroup() |>
      filter(event_id == my_ev) |>
      select(sample_id, PSI),
    by = "sample_id"
  ),
  left_join(
    sim_sf[[my_rep]] |>
      filter(transcript_id == "STRG.6313.3") |>
      select(sample_id, TPM_C = TPM),
    sim_quantifs[[my_rep]] |>
      ungroup() |>
      filter(event_id == my_ev) |>
      select(sample_id, PSI),
    by = "sample_id"
  ),
  by = c("sample_id", "PSI")
) |>
  plotly::plot_ly(x = ~TPM_D, y = ~TPM_C, z = ~PSI)



# Sparse LASSO

# sim1
quantifs_filtered <- sim_quantifs[[1]]
mat_sf_expression <- sim_mat_sf[[1]]

reg_lasso1 <- regression_wrapper(my_ev = my_ev,
                                 regression_method = "lasso",
                                 column = "dPSI_nat",
                                 shuffle = FALSE,
                                 mat_sf_expression = mat_sf_expression,
                                 quants = quantifs_filtered,
                                 intercept = FALSE)

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
true_coefs1 <- sim_true_coefs[[1]] |>
  filter(event_id == my_ev) |>
  mutate(true_coef = if_else(contribution == "inclusion",
                             true_coef,
                             - true_coef)) |>
  group_by(event_id, transcript_id) |>
  summarize(true_coef = sum(true_coef),
            .groups = 'drop')
# mutate(true_coef = if_else(abs(true_coef) > 0.8,
#                            true_coef, 0))

true_coefs1[true_coefs1$true_coef != 0,]



# sim 2
quantifs_filtered <- sim_quantifs[[2]]
mat_sf_expression <- sim_mat_sf[[2]]
reg_lasso2 <- regression_wrapper(my_ev = my_ev, regression_method = "lasso", column = "dPSI_nat",
                                 shuffle = FALSE, mat_sf_expression = mat_sf_expression,
                                 quants = quantifs_filtered,
                                 intercept = FALSE)

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
true_coefs2 <- sim_true_coefs[[2]] |>
  filter(event_id == my_ev) |>
  mutate(true_coef = if_else(contribution == "inclusion",
                             true_coef,
                             - true_coef)) |>
  group_by(event_id, transcript_id) |>
  summarize(true_coef = sum(true_coef),
            .groups = 'drop')
# mutate(true_coef = if_else(abs(true_coef) > 0.8,
#                            true_coef, 0))

true_coefs2[true_coefs2$true_coef != 0,]






comp_to_true <- list(sim1    = reg_lasso1$coefs_sf |> 
                       add_column(sim_rep = "sim1"),
                     sim2 = reg_lasso2$coefs_sf |> 
                       add_column(sim_rep = "sim2")) |>
  list_rbind() |>
  add_column(event_id = my_ev) |>
  left_join(bind_rows(true_coefs1 |> add_column(sim_rep = "sim1"),
                      true_coefs2 |> add_column(sim_rep = "sim2")),
            by = c("event_id", "transcript_id", "sim_rep"))


comp_to_true |>
  filter(transcript_id != "(Intercept)") |>
  ggplot() +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'grey', alpha = 0.5) +
  geom_vline(aes(xintercept = 0), color = 'grey', alpha = 0.5) +
  geom_point(aes(x = true_coef, y = s1)) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'grey90', linetype = 'dashed') +
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
                       ~ do_regression(my_ev = my_ev, regression_method = "lasso", column = "PSI",
                                       shuffle = FALSE,
                                       mat_sf_expression = sim_mat_sf[[.x]],
                                       quants = sim_quantifs[[.x]],
                                       intercept = TRUE) |>
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


comp_to_true |>
  ggplot() +
  theme_classic() +
  geom_hex(aes(x = true_coef, y = s1)) +
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
         FDR = FP/(FP+TP)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = 1-FDR, y = TPR)) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'grey', alpha = .5, linetype = 'dashed')


#~  all replicates, all events ----

get_rsquared <- function(.method, .unit){
  map(sample(events_to_keep, 10),
      \(my_ev){
        reg_lasso_coefs <- map(sample(100, 2),
                               ~ do_regression(my_ev = my_ev, regression_method = .method,
                                               column = .unit,
                                               shuffle = FALSE,
                                               mat_sf_expression = sim_mat_sf[[.x]],
                                               quants = sim_quantifs[[.x]],
                                               intercept = if_else(.unit == "PSI",
                                                                   TRUE, FALSE)) |>
                                 pluck("coefs_sf",
                                       .default = tibble(transcript_id = NA_character_,
                                                         s1 = NA_real_)) |> 
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
  unnest(res) |>
  rename(rsquared = res) |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = regression_method, fill = unit, y = rsquared)) +
  scale_y_continuous(limits = c(0,1))


# this was all using Rsquared as output, what if using the TPR and FDR


get_TPR_FDR <- possibly(function(.method, .unit){
  map(sample(events_to_keep, 10),
      \(my_ev){
        reg_lasso_coefs <- map(sample(100, 10),
                               ~ do_regression(my_ev = my_ev, regression_method = .method,
                                               column = .unit,
                                               shuffle = FALSE,
                                               mat_sf_expression = sim_mat_sf[[.x]],
                                               quants = sim_quantifs[[.x]],
                                               intercept = if_else(.unit == "PSI",
                                                                   TRUE, FALSE)) |>
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
}, otherwise = tibble())


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





