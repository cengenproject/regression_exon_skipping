# Illustration PSI vs TPM for known cases


library(tidyverse)
library(wbData)


tx2g <- wb_load_tx2gene(281)
gids <- wb_load_gene_ids(281) |>
  add_row(X="0", gene_id = "(Intercept)", symbol ="(Intercept)",
          sequence = "(Intercept)", status="Live",biotype="none",name="(Intercept)")




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






## Examples ----


quantifs_filtered <- quantifs_filtered |>
  mutate(target_id = convert_event2_gene_id(event_id),
         target_name = i2s(target_id, gids, TRUE))


# ret-1: hrpr-1,sfa-1,uaf-2,snr-1,prp-38, unc-75

quantifs_filtered |>
  filter(target_name == "ret-1") |> # pull(event_id) |> unique()
  left_join(sf_expression |> filter(gene_name == "hrpr-1") |> group_by(sample_id) |> summarize(TPM = sum(TPM)),
            by = "sample_id") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = TPM, y = PSI, color = neuron_id)) +
  theme(legend.position = "none")

quantifs_filtered |>
  filter(target_name == "ret-1") |> # pull(event_id) |> unique()
  left_join(sf_expression |> filter(gene_name == "sfa-1") |> group_by(sample_id) |> summarize(TPM = sum(TPM)),
            by = "sample_id") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = TPM, y = PSI, color = neuron_id)) +
  theme(legend.position = "none")

quantifs_filtered |>
  filter(target_name == "ret-1") |> # pull(event_id) |> unique()
  left_join(sf_expression |> filter(gene_name == "uaf-2") |> group_by(sample_id) |> summarize(TPM = sum(TPM)),
            by = "sample_id") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = TPM, y = PSI, color = neuron_id)) +
  theme(legend.position = "none")

quantifs_filtered |>
  filter(target_name == "ret-1") |> # pull(event_id) |> unique()
  left_join(sf_expression |> filter(gene_name == "snr-1") |> group_by(sample_id) |> summarize(TPM = sum(TPM)),
            by = "sample_id") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = TPM, y = PSI, color = neuron_id)) +
  theme(legend.position = "none")

quantifs_filtered |>
  filter(target_name == "ret-1") |> # pull(event_id) |> unique()
  left_join(sf_expression |> filter(gene_name == "unc-75") |> group_by(sample_id) |> summarize(TPM = sum(TPM)),
            by = "sample_id") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = TPM, y = PSI, color = neuron_id)) +
  theme(legend.position = "none")


xx <- quantifs_filtered |>
  filter(target_name == "ret-1") |> # pull(event_id) |> unique()
  left_join(sf_expression |> filter(gene_name == "hrpr-1") |> group_by(sample_id) |> summarize(TPM = sum(TPM)),
            by = "sample_id") |>
  mutate(TPM = log(TPM))

plot(xx$TPM, xx$PSI)
lm(PSI ~ TPM, data = xx)
cor.test(xx$TPM, xx$PSI)



xx <- quantifs_filtered |>
  filter(target_name == "ret-1") |> # pull(event_id) |> unique()
  left_join(sf_expression |>
              filter(gene_name %in% c("hrpr-1", "sfa-1", "uaf-2", "snr-1", "prp-38", "unc-75")) |>
              group_by(sample_id, gene_id) |>
              summarize(TPM = sum(TPM)),
            by = "sample_id")
samples <- xx$sample_id |> unique()
y <- quantifs_filtered$PSI[quantifs_filtered$sample_id == samples]
glmnet::cv.glmnet()



# lin-10: hrpr-1

quantifs_filtered |>
  filter(target_name == "lin-10") |> # pull(event_id) |> unique()
  left_join(sf_expression |> filter(gene_name == "hrpr-1") |> group_by(sample_id) |> summarize(TPM = sum(TPM)),
            by = "sample_id") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = TPM, y = PSI, color = neuron_id)) +
  theme(legend.position = "none")


xx <- quantifs_filtered |>
  filter(target_name == "lin-10") |> # pull(event_id) |> unique()
  left_join(sf_expression |> filter(gene_name == "hrpr-1") |> group_by(sample_id) |> summarize(TPM = sum(TPM)),
            by = "sample_id") |>
  mutate(TPM = log(TPM))

mod <- lm(PSI ~ TPM, data = xx)
plot(xx$TPM, xx$PSI); abline(mod)

cor.test(xx$TPM, xx$PSI)



