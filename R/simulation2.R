# Simulation to check whether the permutation test gives results that make sense

# Version 2, expand on the previous version (`simulation_perm_test.R`) for a more realistic simulated dataset

## --- Approach ---
# Keep measured TPMs, just randomize them
# For each event, select 3+3 transcripts non-zero, all other coefficients 0
# For these 6 “true non-0” transcripts, coefficient from Norm(mu=0, sd=1)
# Compute simulation:
#  number of included reads and number of excluded reads simulated from 3 SFs
#  Nincl = Nexcl = Poisson(logistic(sum(coef * TPM)))
#  then PSI = Nincl/(Nincl+Nexcl)
# 
# Create 100 simulated datasets
# 
# Note: in model, use scaled-centered log-TPM



# ----Simulations ----

library(tidyverse)

# # read filtered objects with real data to generate simulations from
# quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
# sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs")

# unfiltered data

quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) |>
  filter(! is.na(PSI))

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") |>
  filter(gene_id %in% putative_splice_factors) |>
  filter(sample_id %in% unique(quantifs$sample_id)) # remove RICr133, Ref, ...



# Parameters of real data ----

counts_measured <- quantifs |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads)) |>
  select(event_id, sample_id, Nincl, Nexcl) |>
  pivot_longer(c(Nincl, Nexcl), names_to = "contribution", values_to = "count") |>
  group_by(event_id, contribution) |>
  nest()

fit <- map(counts_measured$data,
           ~ possibly(\(dat) fitdistrplus::fitdist(data = dat[["count"]],
                                                   distr = "nbinom"))(.x)) |>
  set_names(paste0(counts_measured$event_id,"_",counts_measured$contribution))


all_fits_res <- fit[map_lgl(fit, ~ ! is.null(.x))] |>
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
  facet_wrap(~parameter) +
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
#   mu    meanlog       sdlog
# incl   2.650990   1.430999   
# excl   2.499753   1.637029 

# get whole distribution of size
distrib_size_incl <- all_fits_res |>
  filter(contribution == "Nincl") |>
  pull("size")
distrib_size_excl <- all_fits_res |>
  filter(contribution == "Nexcl") |>
  pull("size")


real_data_fit <- list(mu_incl = fit_mu_incl$estimate,
                      mu_excl = fit_mu_excl$estimate,
                      size_incl = distrib_size_incl,
                      size_excl = distrib_size_excl)
qs::qsave(real_data_fit,
          "data/intermediates/230411_real_data_fit_mu.qs")




# simulate single ----

# prepare simulation data
real_data_fit <- qs::qread("data/intermediates/230411_real_data_fit_mu.qs")

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
  select(transcript_id, sample_id) |>
  mutate(TPM = sample(sf_expression$TPM),
         neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{2,4}$")[,2])

true_coefs <- expand_grid(event_id = unique(quantifs$event_id),
                          contribution = c("inclusion", "exclusion"),
                          transcript_id = unique(sim_sf$transcript_id)) |>
  group_by(event_id, contribution) |>
  nest() |>
  mutate(data = map(data, ~ add_column(.x, true_coef = sample(c(rep(0, nb_tx - 3), rnorm(3)))))) |>
  unnest(data) |> ungroup()

# get list of coefficients for each sample, event, transcript, contribution
randomized_samples <- quantifs |>
  select(event_id, sample_id) |>
  full_join(sim_sf, by = "sample_id", relationship = "many-to-many") |>
  left_join(true_coefs, by = c("event_id", "transcript_id"), relationship = "many-to-many")


# compute

rescale_distr <- function(x, targets){
  targets[["sdlog"]] * x/sd(x) + targets[["meanlog"]]
}

sim_quantifs <- randomized_samples |>
  group_by(event_id, sample_id, contribution) |>
  summarize(val = sum((true_coef*log1p(TPM))),
            .groups = 'drop') |>
  mutate(val_centerd_scaled = if_else(contribution == "inclusion",
                                      rescale_distr(val, real_data_fit$mu_incl),
                                      rescale_distr(val, real_data_fit$mu_excl)),
         val_lnorm = exp(val_centerd_scaled),
         N = rnbinom(n = 2*nb_datapoints,
                     size = if_else(contribution == "inclusion",
                                    sample(real_data_fit$size_incl, size = 2*nb_datapoints, replace = TRUE),
                                    sample(real_data_fit$size_excl, size = 2*nb_datapoints, replace = TRUE)),
                     mu = val_lnorm)) |>
  pivot_wider(id_cols = c(event_id, sample_id),
              names_from = "contribution",
              values_from = "N") |>
  mutate(inclusion = inclusion,
         PSI = inclusion / (inclusion + exclusion),
         nb_reads = inclusion + exclusion,
         neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{2,4}$")[,2])


# Compare histograms with real data
tibble(type = c(rep("measured", nb_datapoints), rep("simul", nb_datapoints)),
       PSI = c(quantifs$PSI, sim_quantifs$PSI)) |>
  ggplot() + theme_classic() +
  geom_density(aes(x = PSI, fill = type), alpha = .3, bw = .01)

tibble(type = c(rep("measured", nb_datapoints), rep("simul", nb_datapoints)),
       PSI = c(quantifs$PSI, sim_quantifs$PSI)) |>
  ggplot() + theme_classic() +
  geom_freqpoly(aes(x = PSI, color = type), bins = 100)


tibble(type = c(rep("measured", nb_datapoints), rep("simul", nb_datapoints)),
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

sample2neuron <- function(sample_ids){
  res <- str_extract(sample_ids, "^([A-Z0-9ef]{2,4})r\\d{1,3}$", group = 1L)
  if(any(is.na(res))) warning(sum(is.na(res)), " NAs returned: ", head(sample_ids[is.na(res)]))
  res
}

# Filter neurons ----

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



# Re-Make tx matrix

mat_tx_expr_real <- sf_expression |>
  pivot_wider(id_cols = transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr_real <- mat_tx_expr_real[rowSums(mat_tx_expr_real) != 0,]
dim(mat_tx_expr_real)


mat_tx_expr_sim <- sim_sf |>
  pivot_wider(id_cols = transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr_sim <- mat_tx_expr_sim[rowSums(mat_tx_expr_sim) != 0,]
dim(mat_tx_expr_sim)





#~ Filter samples based on tx data ----

# PCA
pc_real <- prcomp(t(log1p(mat_tx_expr_real)), center = TRUE, scale. = TRUE)

pcdf <- cbind(as.data.frame(pc_real$x), data.frame(sample_id = colnames(mat_tx_expr_real))) |>
  mutate(neuron_id = sample2neuron(sample_id))

plot(pcdf[,1:2])
ggplot(pcdf, aes(x = PC1 , y = PC2, color = neuron_id)) +
  theme_classic() +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id)) +
  geom_abline(aes(intercept = -150, slope = -.3)) +
  hues::scale_color_iwanthue()

keep_samples_real <- unique(sf_expression$sample_id) |>
  setdiff(c("RICr133","PVMr122","AVKr113"))


pc_sim <- prcomp(t(log1p(mat_tx_expr_sim)), center = TRUE, scale. = TRUE)

pcdf <- cbind(as.data.frame(pc_sim$x), data.frame(sample_id = colnames(mat_tx_expr_sim))) |>
  mutate(neuron_id = sample2neuron(sample_id))

plot(pcdf[,1:2])
ggplot(pcdf, aes(x = PC1 , y = PC2, color = neuron_id)) +
  theme_classic() +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id)) +
  geom_abline(aes(intercept = -150, slope = -.3)) +
  hues::scale_color_iwanthue()

keep_samples_sim <- unique(sim_sf$sample_id) |>
  setdiff(c("RICr133","PVMr122","AVKr113"))





# Filter bad samples
sf_expression <- sf_expression |>
  filter(sample_id %in% keep_samples_real)

quantifs <- quantifs |>
  filter(sample_id %in% keep_samples_real)

sim_sf <- sim_sf |>
  filter(sample_id %in% keep_samples_sim)

sim_quantifs <- sim_quantifs |>
  filter(sample_id %in% keep_samples_sim)




# Re-Make tx matrix
mat_tx_expr_real <- sf_expression |>
  pivot_wider(id_cols = transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr_real <- mat_tx_expr_real[rowSums(mat_tx_expr_real) != 0,]
dim(mat_tx_expr_real)


mat_tx_expr_sim <- sim_sf |>
  pivot_wider(id_cols = transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr_sim <- mat_tx_expr_sim[rowSums(mat_tx_expr_sim) != 0,]
dim(mat_tx_expr_sim)







#~ Filter samples based on exon data ----

# make matrix
mat_quantifs_real <- quantifs |>
  pivot_wider(id_cols = event_id,
              names_from = "sample_id",
              values_from = "PSI") |>
  column_to_rownames("event_id") |>
  as.matrix()
dim(mat_quantifs_real)

mat_quantifs_sim <- sim_quantifs |>
  pivot_wider(id_cols = event_id,
              names_from = "sample_id",
              values_from = "PSI") |>
  column_to_rownames("event_id") |>
  as.matrix()
dim(mat_quantifs_sim)



# need to tackle NAs before PCA
hist(rowSums(is.na(mat_quantifs_real)))
hist(colSums(is.na(mat_quantifs_real)))

mat_quantifs_real <- mat_quantifs_real[rowSums(is.na(mat_quantifs_real)) < 100,
                             colSums(is.na(mat_quantifs_real)) < 500]
dim(mat_quantifs_real)

mat_quantifs_real <- missMDA::imputePCA(mat_quantifs_real)
dim(mat_quantifs_real$completeObs)



hist(rowSums(is.na(mat_quantifs_sim)))
hist(colSums(is.na(mat_quantifs_sim)))

mat_quantifs_sim <- mat_quantifs_sim[rowSums(is.na(mat_quantifs_sim)) < 100,
                             colSums(is.na(mat_quantifs_sim)) < 500]
dim(mat_quantifs_sim)

mat_quantifs_sim <- missMDA::imputePCA(mat_quantifs_sim)
dim(mat_quantifs_sim$completeObs)


# PCA
pc_real <- prcomp(t(mat_quantifs_real$completeObs), center = TRUE, scale. = TRUE)
pcdf <- cbind(as.data.frame(pc_real$x), data.frame(sample_id = colnames(mat_quantifs_real$completeObs))) |>
  mutate(neuron_id = sample2neuron(sample_id))

plot(pcdf[,1:2])
ggplot(pcdf, aes(x = PC1 , y = PC2, color = neuron_id)) +
  theme_classic() +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id)) +
  hues::scale_color_iwanthue()


pc_sim <- prcomp(t(mat_quantifs_sim$completeObs), center = TRUE, scale. = TRUE)
pcdf <- cbind(as.data.frame(pc_sim$x), data.frame(sample_id = colnames(mat_quantifs_sim$completeObs))) |>
  mutate(neuron_id = sample2neuron(sample_id))

plot(pcdf[,1:2])
ggplot(pcdf, aes(x = PC1 , y = PC2, color = neuron_id)) +
  theme_classic() +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id)) +
  hues::scale_color_iwanthue()






# Filter tx ----

# PCA on tx
pc_real <- prcomp(log1p(mat_tx_expr_real), center = TRUE, scale. = TRUE)

pcdf <- cbind(as.data.frame(pc_real$x), data.frame(tx_id = rownames(mat_tx_expr_real)))
plot(pcdf[,1:2])
abline(v = 80)
pcdf[pcdf$PC1 > 60,"tx_id"] |> writeClipboard()
#> ribosomal, mitochondrial, ncRNAs, etc

# table(wbData::wb_tx2g(pcdf[pcdf$PC1 > 60,"tx_id"], tx2g) %in% wormDatasets::worm_putative_splice_factors$gene_id)
keep_tx_real <- rownames(mat_tx_expr_real) |>
  setdiff(pcdf[pcdf$PC1 > 60,"tx_id"])





pc_sim <- prcomp(log1p(mat_tx_expr_sim), center = TRUE, scale. = TRUE)

pcdf <- cbind(as.data.frame(pc_sim$x), data.frame(tx_id = rownames(mat_tx_expr_sim)))
plot(pcdf[,1:2])
abline(v = 80)
pcdf[pcdf$PC1 > 60,"tx_id"] |> writeClipboard()
#> ribosomal, mitochondrial, ncRNAs, etc


# table(wbData::wb_tx2g(pcdf[pcdf$PC1 > 60,"tx_id"], tx2g) %in% wormDatasets::worm_putative_splice_factors$gene_id)

keep_tx_sim <- rownames(mat_tx_expr_sim) |>
  setdiff(pcdf[pcdf$PC1 > 60,"tx_id"])



# Filter bad tx
sf_expression <- sf_expression |>
  filter(transcript_id %in% keep_tx_real)
sim_sf <- sim_sf |>
  filter(transcript_id %in% keep_tx_sim)



# Re-Make tx matrix
mat_tx_expr_real <- sf_expression |>
  pivot_wider(id_cols = transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr_real <- mat_tx_expr_real[rowSums(mat_tx_expr_real) != 0,]
dim(mat_tx_expr_real)


mat_tx_expr_sim <- sim_sf |>
  pivot_wider(id_cols = transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr_sim <- mat_tx_expr_sim[rowSums(mat_tx_expr_sim) != 0,]
dim(mat_tx_expr_sim)






#~ Filter exons ----

# make matrix
mat_quantifs_real <- quantifs |>
  pivot_wider(id_cols = event_id,
              names_from = "sample_id",
              values_from = "PSI") |>
  column_to_rownames("event_id") |>
  as.matrix()
dim(mat_quantifs_real)


mat_quantifs_sim <- sim_quantifs |>
  pivot_wider(id_cols = event_id,
              names_from = "sample_id",
              values_from = "PSI") |>
  column_to_rownames("event_id") |>
  as.matrix()
dim(mat_quantifs_sim)



# need to tackle NAs before PCA
hist(rowSums(is.na(mat_quantifs_real)))
hist(colSums(is.na(mat_quantifs_real)))

mat_quantifs_real <- mat_quantifs_real[rowSums(is.na(mat_quantifs_real)) < 100,
                             colSums(is.na(mat_quantifs_real)) < 500]
dim(mat_quantifs_real)

mat_quantifs_real <- missMDA::imputePCA(mat_quantifs_real)
dim(mat_quantifs_real$completeObs)




hist(rowSums(is.na(mat_quantifs_sim)))
hist(colSums(is.na(mat_quantifs_sim)))

mat_quantifs_sim <- mat_quantifs_sim[rowSums(is.na(mat_quantifs_sim)) < 100,
                             colSums(is.na(mat_quantifs_sim)) < 500]
dim(mat_quantifs_sim)

mat_quantifs_sim <- missMDA::imputePCA(mat_quantifs_sim)
dim(mat_quantifs_sim$completeObs)



# PCA
pc_real <- prcomp(mat_quantifs_real$completeObs, center = TRUE, scale. = TRUE)

pcdf <- cbind(as.data.frame(pc_real$x), data.frame(sample_id = rownames(mat_quantifs_real$completeObs)))
plot(pcdf[,1:2])



pc_sim <- prcomp(mat_quantifs_sim$completeObs, center = TRUE, scale. = TRUE)

pcdf <- cbind(as.data.frame(pc_sim$x), data.frame(sample_id = rownames(mat_quantifs_sim$completeObs)))
plot(pcdf[,1:2])





# Filter events ----

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
  geom_hline(aes(yintercept = 0.05), color = 'grey')

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


events_to_keep_real <- intersect(events_to_keep_n_samples_real,
                            events_to_keep_variability_real)

quantifs_filtered_real <- quantifs_filtered_n_reads_real |>
  filter(event_id %in% events_to_keep_real) |>
  filter(! is.na(PSI))

all.equal(sort(events_to_keep_real), sort(unique(quantifs_filtered_real$event_id)))


quantifs_filtered_real |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads),color='grey', bins = 100) +
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
  geom_hline(aes(yintercept = 0.05), color = 'grey')

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




events_to_keep_sim <- intersect(events_to_keep_n_samples_sim,
                            events_to_keep_variability_sim)

quantifs_filtered_sim <- quantifs_filtered_n_reads_sim |>
  filter(event_id %in% events_to_keep_sim) |>
  filter(! is.na(PSI))

all.equal(sort(events_to_keep_sim), sort(unique(quantifs_filtered_sim$event_id)))


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







# qs::qsave(quantifs_filtered_sim, "data/intermediates/230411_simulation/quantifs_filtered.qs")
# qs::qsave(sim_sf, "data/intermediates/230411_simulation/sim_sf.qs")
# qs::qsave(true_coefs, "data/intermediates/230411_simulation/true_coefs.qs")





# Check simulated PSI
(my_ev <- sample(events_to_keep_sim, 1))
bind_rows(
  quantifs_filtered_real |>
    filter(event_id == my_ev) |>
    select(ends_with("_id"), PSI) |>
    add_column(type = "measured"),
  quantifs_filtered_sim |>
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

hist(quantifs_filtered_real$PSI, breaks = 150)
hist(quantifs_filtered_sim$PSI, breaks = 150)

bind_rows(
  quantifs_filtered_real |>
    select(ends_with("_id"), PSI) |>
    add_column(type = "measured"),
  quantifs_filtered_sim |>
    select(ends_with("_id"), PSI) |>
    mutate(neuron_id = str_match(sample_id, "^([A-Z0-9]+)r[0-9]+$")[,2]) |>
    add_column(type = "simulated")) |>
  ggplot() +
  theme_classic() +
  geom_freqpoly(aes(x = PSI, color = type), bins = 100)










# Simulate 100 datasets ----


library(tidyverse)

# unfiltered data
quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) |>
  filter(! is.na(PSI))

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") |>
  filter(gene_id %in% putative_splice_factors) |>
  filter(sample_id %in% unique(quantifs$sample_id)) # remove RICr133, Ref, ...

real_data_fit <- qs::qread("data/intermediates/230411_real_data_fit_mu.qs")




rescale_distr <- function(x, targets){
  targets[["sdlog"]] * x/sd(x) + targets[["meanlog"]]
}

simulate_single <- function(real_data_fit, sf_expression, quantifs){
  nb_tx  <- sf_expression |>
    pull(transcript_id) |>
    unique() |>
    length()
  
  nb_datapoints <- quantifs |> nrow()
  
  sim_sf <- sf_expression |>
    select(transcript_id, sample_id) |>
    mutate(TPM = sample(sf_expression$TPM),
           neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{2,4}$")[,2])
  
  true_coefs <- expand_grid(event_id = unique(quantifs$event_id),
                            contribution = c("inclusion", "exclusion"),
                            transcript_id = unique(sim_sf$transcript_id)) |>
    group_by(event_id, contribution) |>
    nest() |>
    mutate(data = map(data, ~ add_column(.x, true_coef = sample(c(rep(0, nb_tx - 3), rnorm(3)))))) |>
    unnest(data) |> ungroup()
  
  # get list of coefficients for each sample, event, transcript, contribution
  randomized_samples <- quantifs |>
    select(event_id, sample_id) |>
    full_join(sim_sf, by = "sample_id", relationship = "many-to-many") |>
    left_join(true_coefs, by = c("event_id", "transcript_id"), relationship = "many-to-many")
  
  
  # compute
  sim_quantifs <- randomized_samples |>
    group_by(event_id, sample_id, contribution) |>
    summarize(val = sum((true_coef*log1p(TPM))),
              .groups = 'drop') |>
    mutate(val_centerd_scaled = if_else(contribution == "inclusion",
                                        rescale_distr(val, real_data_fit$mu_incl),
                                        rescale_distr(val, real_data_fit$mu_excl)),
           val_lnorm = exp(val_centerd_scaled),
           N = rnbinom(n = 2*nb_datapoints,
                       size = if_else(contribution == "inclusion",
                                      sample(real_data_fit$size_incl, size = 2*nb_datapoints, replace = TRUE),
                                      sample(real_data_fit$size_excl, size = 2*nb_datapoints, replace = TRUE)),
                       mu = val_lnorm)) |>
    pivot_wider(id_cols = c(event_id, sample_id),
                names_from = "contribution",
                values_from = "N") |>
    mutate(inclusion = inclusion,
           PSI = inclusion / (inclusion + exclusion),
           nb_reads = inclusion + exclusion,
           neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{2,4}$")[,2])
  
  list(sim_quantifs = sim_quantifs, sim_sf = sim_sf, true_coefs = true_coefs)
}

# 20s per replicate, 30' for 100
sim_replicated <- replicate(100,
                            simulate_single(real_data_fit, sf_expression, quantifs),
                            simplify = FALSE)


# qs::qsave(sim_replicated, "data/intermediates/230411_simulation/rep_simulations.qs")







# ---- Regression ----

library(tidyverse)
library(glmnet)

source("R/regression_functions.R")


# Read data ----
sim_replicated <- qs::qread("data/intermediates/230411_simulation/sim_sf.qs")

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





