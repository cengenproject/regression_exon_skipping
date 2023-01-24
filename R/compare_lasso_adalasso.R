# Compare LASSO and adaptive LASSO (see scripts "explore_[adaptive]_lasso_on_PSI_vs_TPM")


library(tidyverse)
library(glmnet)
library(wbData)


tx2g <- wb_load_tx2gene(281)
gids <- wb_load_gene_ids(281)


# Data ----

events_coordinates <- read_tsv("data/export_for_arman/221111_events_coordinates.tsv")

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)


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




sf_targets <- jsonlite::read_json("data/export_for_arman/sf_targets_v2.json")

sf2target <- tibble(sf_id = map_chr(sf_targets, \(x) x[["SF"]]),
                    target_id = map(sf_targets, \(x) x[["targets"]])) |>
  mutate(sf_id = if_else(sf_id == "mec-8 ad", "mec-8", sf_id),
         sf_id = map_chr(sf_id,
                         \(x) `if`(startsWith(x, "WBGene"), x, s2i(x, gids, warn_missing = TRUE)))) |>
  unnest_longer(target_id) |>
  filter(target_id != "0") |>
  mutate(sf_name = i2s(sf_id, gids, warn_missing = FALSE),
         target_name = i2s(target_id, gids, warn_missing = TRUE)) |>
  distinct()





# Load cached result from other scripts
reg_lasso <- qs::qread("data/intermediates/230117_replicated_regression_log_cache_filt.qs")
reg_adalasso <- qs::qread("data/intermediates/230118_replicated_regression_adalasso_cache.qs")

events_names <- readLines("data/intermediates/230119_events_to_keep.txt")


# Extract coefficient values
coefs_lasso <- imap(reg_lasso,
                    \(replicate, ind) {
                      map2(replicate, events_names,
                           \(rep, .event_id) {
                             pluck(rep, "coefs_sf", .default = tibble()) |>
                               add_column(event_id = .event_id)
                           }) |>
                        list_rbind() |>
                        add_column(replicate = paste0("replicate_", ind))
                    }) |>
  list_rbind()

coefs_adalasso <- imap(reg_adalasso,
                    \(replicate, ind) {
                      map2(replicate, events_names,
                           \(rep, .event_id) {
                             pluck(rep, "coefs_sf", .default = tibble()) |>
                               add_column(event_id = .event_id)
                           }) |>
                        list_rbind() |>
                        add_column(replicate = paste0("replicate_", ind))
                    }) |>
  list_rbind()


# Find stable coefficients ----
coefs_stability <- bind_rows(coefs_lasso |> add_column(method = "lasso"),
                             coefs_adalasso |> add_column(method = "adalasso")) |>
  group_by(method, event_id, transcript_id) |>
  summarize(median = median(s1),
            mean = mean(s1),
            mad = mad(s1),
            sd = sd(s1),
            .groups = "drop")

coefs_stability <- coefs_stability |>
  filter(transcript_id != "(Intercept)")


table(coefs_stability$mean == 0, coefs_stability$median == 0)

coefs_stability_nonnull <- coefs_stability |>
  filter(median != 0)


# qs::qsave(coefs_stability_nonnull, "data/intermediates/230119_coefs_stability_nonnull_cache.qs")


coefs_stability_nonnull <- qs::qread("data/intermediates/230119_coefs_stability_nonnull_cache.qs")


# Filter stable events

ggplot(coefs_stability_nonnull) +
  theme_classic() +
  geom_point(aes(x = abs(median), y = mad), alpha = .2) +
  scale_x_log10() + scale_y_log10() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_wrap(~method)

ggplot(coefs_stability_nonnull) +
  theme_classic() +
  geom_point(aes(x = abs(mean), y = sd), alpha = .2) +
  scale_x_log10() + scale_y_log10() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_wrap(~method)





# add target name, add whether known interaction
coefs_stability_withknown <- coefs_stability_nonnull |>
  left_join(events_coordinates |>
              select(event_id, gene_id),
            by = "event_id") |>
  rename(target_id = gene_id) |>
  mutate(sf_id = convert_sf_tx2g(transcript_id),
         is_computed = TRUE) |>
  # add as lasso and adalasso separately
  full_join(sf2target |>
              select(sf_id, target_id) |>
              add_column(is_known = TRUE,
                         method = "adalasso"),
            by = c("method", "target_id", "sf_id")) |>
  full_join(sf2target |>
              select(sf_id, target_id) |>
              add_column(is_known = TRUE,
                         method = "lasso"),
            by = c("method", "target_id", "sf_id")) |>
  mutate(is_known = is_known.x | is_known.y,
         is_known = !is.na(is_known),
         is_computed = !is.na(is_computed)) |>
  select(method,
         target_event_id = event_id, target_gene_id = target_id,
         sf_gene_id = sf_id, sf_tx_id = transcript_id,
         is_computed, is_known,
         mean_computed_interaction = mean,
         sd_computed_interaction = sd)


table(coefs_stability_withknown$method, coefs_stability_withknown$is_computed,
      useNA = 'ifany', dnn=c("method", "computed")) |> addmargins()


# try a series of thresholds

coefs_stability_withknown |>
  filter(is_computed) |>
  group_by(method, is_known) |>
  mutate(thresholds = list((1:100)/1000)) |>
  select(method, is_known, mean_computed_interaction, thresholds) |>
  unnest(c(thresholds)) |>
  mutate(kept = abs(mean_computed_interaction) >= thresholds) |>
  group_by(thresholds, .add = TRUE) |>
  summarize(nb_kept = sum(kept)) |>
  mutate(nb_kept = round(100*nb_kept/max(nb_kept))) |>
  ungroup() |>
  mutate(is_known = if_else(is_known, "known", "new")) |>
  pivot_wider(c(method, thresholds),
              names_from = "is_known",
              values_from = "nb_kept") |>
  ggplot() +
  theme_classic() +
  geom_line(aes(x = new, y = known, color = method)) +
  geom_point(aes(x = new, y = known, color = method)) +
  ggrepel::geom_text_repel(aes(x = new, y = known, color = method, label = thresholds)) +
  geom_vline(aes(xintercept = 26), color = 'grey', linetype = 'dashed')



coefs_stability_withknown |>
  filter(is_computed) |>
  group_by(method, is_known) |>
  mutate(thresholds = list((1:100)/1000)) |>
  select(method, is_known, mean_computed_interaction, thresholds) |>
  unnest(c(thresholds)) |>
  mutate(kept = abs(mean_computed_interaction) >= thresholds) |>
  group_by(thresholds, .add = TRUE) |>
  summarize(nb_kept = sum(kept)) |>
  mutate(nb_kept_scaled = nb_kept/max(nb_kept)) |>
  ungroup() |>
  ggplot() +
  theme_classic() +
  geom_line(aes(x = thresholds, y = nb_kept_scaled, color = is_known)) +
  facet_wrap(~ method) +
  geom_vline(aes(xintercept = 0.015), color = 'grey', linetype = 'dashed')


coefs_stability_withknown |>
  filter(is_computed) |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = abs(mean_computed_interaction), color = is_known), bw = .07) +
  scale_x_log10() +
  # scale_x_continuous(limits=c(.002,.02)) +
  facet_wrap(~method) +
  geom_vline(aes(xintercept = .015), linetype = 'dotted')





coefs_stability_withknown <- coefs_stability_withknown |>
  mutate(is_stable = abs(mean_computed_interaction) >= 0.015 &
           abs(mean_computed_interaction) >= 0.5 * sd_computed_interaction,
         is_stable = if_else(!is.na(is_stable), is_stable, FALSE))



coefs_stability_withknown |>
  filter(is_computed) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = abs(mean_computed_interaction), y = sd_computed_interaction,
                 color = is_known, shape = is_stable),
             alpha = .1) +
  scale_x_log10() + scale_y_log10() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_wrap(~method) +
  # geom_abline(aes(intercept = -log10(.02), slope = 2), linetype = 'dashed') +
  geom_vline(aes(xintercept = .015), linetype = 'dotted') +
  scale_color_manual(values = c("black","red"))



my_line <- function(x) 2*x-0.015

coefs_stability_withknown |>
  filter(is_computed) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = abs(mean_computed_interaction), y = sd_computed_interaction,
                 color = is_stable),
             alpha = .2) +
  scale_x_log10() + scale_y_log10() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_wrap(~method) +
  geom_function(fun = my_line) +
  # geom_abline(aes(intercept = -log10(.015), slope = 2), linetype = 'dashed') +
  geom_vline(aes(xintercept = .015), linetype = 'dotted') +
  scale_color_manual(values = c("black", "red"))



#~ Compare with known ----

# collapse several events in same target gene or several tx from same SF gene
coefs_collapsed <- coefs_stability_withknown |>
  select(-mean_computed_interaction, -sd_computed_interaction) |>
  pivot_wider(names_from = "method",
              values_from = "is_stable") |>
  group_by(target_gene_id, sf_gene_id,) |>
  summarize(is_known = any(is_known),
            adalasso = any(adalasso),
            lasso = any(lasso),
            .groups = 'drop')


coefs_collapsed |>
  summarize(nb_known = sum(is_known),
            nb_tested_adalasso   = sum(!is.na(adalasso)),
            nb_stable_adalasso   = sum(adalasso, na.rm = TRUE),
            nb_overlaps_adalasso = sum(is_known & adalasso, na.rm = TRUE),
            nb_tested_lasso = sum(!is.na(lasso)),
            nb_stable_lasso = sum(lasso, na.rm = TRUE),
            nb_overlaps_lasso = sum(is_known & lasso, na.rm = TRUE)) |>
  as.data.frame()

# AdaLASSO: 2009 computed as stable, 316+4836 known but not found, 137 overlap




intersected_stable_coefs <- coefs_stability_withknown |>
  filter(is_stable) |>
  rename(computed_sf_tx_id = transcript_id) |>
  # nest computed sf list
  mutate(computed_sf_gene_id = convert_sf_tx2g(computed_sf_tx_id)) |>
  group_by(method, event_id, target_id) |>
  nest(computed = c(computed_sf_tx_id, computed_sf_gene_id)) |>
  ungroup() |>
  group_by(method, event_id, target_id) |>
  nest(known = c(sf_id, sf_name))

# find intersection size of computed and known SF
overlap_intersected_stable_coefs <- intersected_stable_coefs |>
  mutate(nb_sf_computed = map_int(computed, ~length(.x$computed_sf_gene_id)),
         nb_sf_known = map_int(known, ~length(.x$sf_id)),
         nb_sf_overlap = map2_int(intersected_stable_coefs$computed, intersected_stable_coefs$known,
                                      ~ length(intersect(.x$computed_sf_gene_id, .y$sf_id)))
         )


sum(overlap_intersected_stable_coefs$nb_sf_computed)
sum(overlap_intersected_stable_coefs$nb_sf_overlap)
sum(overlap_intersected_stable_coefs$nb_sf_known)
25/482
#> 5% no filtering
178/1768
#> 10% after filtering
274/1768
#> 15% with adaptive lasso


rep_overlap <- pbapply::pbreplicate(500,
                                    sum(map2_int(sample(intersected_coefs$computed), sample(intersected_coefs$known),
                                                 ~ length(intersect(.x$computed_sf_gene_id, .y$sf_id))))
)

hist(rep_overlap, breaks = 30,
     xlab = "Number of overlapping interactions under randomization",
     main = NULL); abline(v = sum(overlap_intersected_coefs$nb_sf_overlap), col = 'red')
table(rep_overlap >= sum(overlap_intersected_coefs$nb_sf_overlap))






