# Created in simulation v9 (probably the same as used for simulation v6-8)
# Load "raw" data as exported from `quantifs_exon_skipping`
# prefilter to keep only variable events etc
# then will be used as basis for simulation

# Note: this is inspired by, and somewhat equivalent to `filter_PSI_and_TPM.R` (but with additional filtering for variable events)
# However, the results of `filter_PSI_and_TPM.R` are NOT used here

#re-run (minor adaptations) 240308 for bsn12 new dataset
# 240429: use all RBPs in addition to putative SFs

# Prefilter ----

library(tidyverse)

quantifs <- read_tsv("../quantif_exon_skipping/data/export_for_arman/240308_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) |>
  filter(! is.na(PSI)) |>
  filter(startsWith(event_id, "SE_"))


eulerr::euler(list(
  putative_sf = wormDatasets::worm_putative_splice_factors |> filter(keep == 1 | keep == 2) |> pull(gene_id),
  tamburino = wormDatasets::list_rbp_tamburino2013$gene_id,
  matiagonzalez = wormDatasets::list_rbp_matiagonzalez2015 |> filter(selected == 1) |> pull(gene_id)
)) |>
  plot()

putative_splice_factors <- union(
  wormDatasets::worm_putative_splice_factors |> filter(keep == 1 | keep == 2) |> pull(gene_id),
  wormDatasets::list_rbp_tamburino2013$gene_id) |>
  union(wormDatasets::list_rbp_matiagonzalez2015 |> filter(selected == 1) |> pull(gene_id))



sf_expression <- read_tsv("data/export_for_arman/231208_t_exp.tsv") |>
  filter(gene_id %in% putative_splice_factors) |>
  filter(sample_id %in% unique(quantifs$sample_id)) # remove RICr133, Ref, ...


# filter SF expression ----

# load Alec's sc-bulk integration
gene_expression <- read.csv("../quantif_exon_skipping/data/thresholded_gene_expression/bsn12_subtracted_integrated_binarized_expression_withVDDD_FDR0.05_030424.tsv",
                            sep = "\t") |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id,
               names_to = "neuron_id",
               values_to = "expression")

sf_expression |>
  left_join(gene_expression,
            by = c("gene_id", "neuron_id")) |>
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = factor(expression), y = log10(TPM + 1)))

noneur <- unique(sf_expression$gene_id)[!unique(sf_expression$gene_id) %in% wormDatasets::genes_by_pattern$present_in_neurons]

sf_expression |>
  mutate(nonneur = gene_id %in% noneur) |>
  ggplot() +theme_classic() +
  geom_density(aes(x = log10(TPM +1), fill = nonneur), alpha = .2)

sf_expression |>
  mutate(nonneur = gene_id %in% noneur) |>
  filter(nonneur & log10(TPM+1)>2) |>
  pull(gene_id) |> unique()
  

sf_expression |>
  mutate(logTPM = log10(TPM+1)) |>
  summarize(mean_TPM = mean(logTPM),
            sd_TPM = sd(logTPM),
            .by = c(transcript_id, gene_id)) |>
  ggplot() +
  theme_classic()+
  geom_point(aes(x = mean_TPM, y = sd_TPM, color = sd_TPM >= .15), alpha = .1) +
  geom_hline(aes(yintercept = .15))


tx_to_keep <- sf_expression |>
  mutate(logTPM = log10(TPM+1)) |>
  summarize(mean_TPM = mean(logTPM),
            sd_TPM = sd(logTPM),
            .by = c(transcript_id, gene_id)) |>
  filter(sd_TPM >= .15) |>
  pull(transcript_id)

# my_gene <- sample(tx_to_keep$gene_id, 1)
# tx_to_keep |> filter(gene_id == my_gene)
# sf_expression |> filter(gene_id == my_gene) |>
#   ggplot() +
#   theme_classic() +
#   ggbeeswarm::geom_quasirandom(aes(x = transcript_id, y = log10(TPM+1)))

sf_expression <- sf_expression |>
  filter(transcript_id %in% tx_to_keep)


#~ Filter neurons ----

# Filter neurons with too few samples
keep_neurons <- sf_expression |>
  select(sample_id, neuron_id) |>
  distinct() |>
  dplyr::count(neuron_id) |>
  filter(n > 2) |>
  pull(neuron_id)


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
(
  quantifs_filtered_n_reads |>
    filter(! is.na(PSI)) |>
    group_by(event_id) |>
    summarize(nb_samples = n(),
              nb_neurons = n_distinct(neuron_id)) |>
    ggplot() + theme_classic() +
    geom_point(aes(x = nb_neurons, y = nb_samples), alpha = .2) +
    xlab("Number of neurons") + ylab("Number of samples") +
    geom_hline(aes(yintercept = 70.5), color = 'darkred') +
    geom_vline(aes(xintercept = 23.5), color = 'darkred')
) |>
  ggExtra::ggMarginal(type = "histogram")



events_to_keep_n_samples <- quantifs_filtered_n_reads |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  filter(nb_samples > 70,
         nb_neurons > 23) |>
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
  geom_hline(aes(yintercept = 0.03), color = 'grey') +
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
  filter(sd_PSI_btw_neurs > 0.03) |>
  pull(event_id)


quantifs_filtered <- quantifs_filtered_nsamples |>
  filter(event_id %in% events_to_keep_variability)


quantifs_filtered |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = nb_reads), color='grey', bins = 100) +
  scale_x_log10()


#~ Save prefiltering ----

qs::qsave(quantifs_filtered,
          "data/graph_power4/inputs/240429_preprocessed_quantifs_filtered.qs")

qs::qsave(sf_expression,
          "data/graph_power4/inputs/240429_preprocessed_sf_expression.qs")



