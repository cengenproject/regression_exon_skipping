# Created in simulation v9 (probably the same as used for simulation v6-8)
# Load "raw" data as exported from `quantifs_exon_skipping`
# prefilter to keep only variable events etc
# then will be used as basis for simulation

# Note: this is inspired by, and somewhat equivalent to `filter_PSI_and_TPM.R` (but with additional filtering for variable events)
# However, the results of `filter_PSI_and_TPM.R` are NOT used here



# Prefilter ----

library(tidyverse)

quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) |>
  filter(! is.na(PSI))

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") |>
  filter(gene_id %in% putative_splice_factors) |>
  filter(sample_id %in% unique(quantifs$sample_id)) # remove RICr133, Ref, ...




#~ Filter neurons ----

# Filter neurons with too few samples. Also remove Ref
keep_neurons <- sf_expression |>
  select(sample_id, neuron_id) |>
  distinct() |>
  dplyr::count(neuron_id) |>
  filter(n > 2) |>
  pull(neuron_id) |>
  setdiff("Ref")


sf_expression <- sf_expression |>
  filter(neuron_id %in% keep_neurons)


quantifs <- quantifs |>
  filter(neuron_id %in% keep_neurons)


#~ Filter events ----

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
  geom_hline(aes(yintercept = 100.5), color = 'darkred') +
  geom_vline(aes(xintercept = 32.5), color = 'darkred')



events_to_keep_n_samples <- quantifs_filtered_n_reads |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  filter(nb_samples > 100,
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
  geom_hline(aes(yintercept = 0.05), color = 'grey') +
  ggtitle("Real data")

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


quantifs_filtered <- quantifs_filtered_nsamples |>
  filter(event_id %in% events_to_keep_variability) |>
  filter(! is.na(PSI))

quantifs_filtered |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads), color='grey', bins = 100) +
  scale_x_log10()


#~ Save prefiltering ----

write_lines(unique(quantifs_filtered$event_id),
            "data/intermediates/230512_simulation_v10/events_to_keep.txt")
write_lines(unique(quantifs_filtered$neuron_id),
            "data/intermediates/230512_simulation_v10/neurons_to_keep.txt")

qs::qsave(quantifs_filtered,
          "data/intermediates/230512_simulation_v10//quantifs_prefiltered.qs")

qs::qsave(sf_expression,
          "data/intermediates/230512_simulation_v10//sf_prefiltered.qs")



