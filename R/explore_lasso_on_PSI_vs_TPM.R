

library(tidyverse)
library(glmnet)
library(wbData)


tx2g <- wb_load_tx2gene(281)
gids <- wb_load_gene_ids(281)
i2s
# Read data ----

quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) 

putative_splice_factors <- wormDatasets::worm_putative_splice_factors |>
  filter(keep == 1 | keep == 2) |>
  pull(gene_id)

sf_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") |>
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



# First pass Sparse Regression ----

sparse_regression <- function(my_ev){
  x <- sf_expression |>
    select(transcript_id, sample_id, TPM) |>
    mutate(TPM = log(TPM + 1)) |>
    pivot_wider(id_cols = sample_id,
                names_from = "transcript_id",
                values_from = "TPM") |>
    column_to_rownames("sample_id") |>
    as.matrix() |>
    scale()
  x <- x[,! apply(x, 2, \(col) any(is.na(col)))]
  
  y <- quantifs[quantifs$event_id == my_ev, c("sample_id", "PSI")] |>
    column_to_rownames("sample_id") |>
    filter(!is.na(PSI)) |>
    as.matrix()
  
  table(rownames(y) %in% rownames(x))
  table(rownames(x) %in% rownames(y))
  
  x <- x[rownames(y),]
  
  
  n <- nrow(x)
  train <- sample(n, round(.7*n))
  
  
  
  cvfit <- cv.glmnet(x[train,], y[train], nfolds = 30)
  
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


# first_pass_regression <- quantifs |>
#   filter(! is.na(PSI)) |>
#   group_by(event_id) |>
#   summarize(nb_samples = n(),
#             nb_neurons = n_distinct(neuron_id)) |>
#   arrange(desc(nb_samples), desc(nb_neurons)) |>
#   mutate(fit = map(event_id, possibly(sparse_regression, otherwise = list(rsquare=NA_real_,nb_coefs=NA_real_)),
#                    .progress = TRUE)) |>
#   mutate(rsquare = map_dbl(fit, \(x) x[["rsquare"]]),
#          nb_coefs = map_dbl(fit, \(x) x[["nb_coefs"]]),
#          prediction_on_test = map(fit, \(x) x[["prediction_on_test"]]),
#          coefs_sf = map(fit, \(x) x[["coefs_sf"]])) |>
#   select(-fit)
# qs::qsave(first_pass_regression, "data/intermediates/230113_fits_for_quantifs_scaled_log_cache.qs")


first_pass_regression <- qs::qread("data/intermediates/230113_fits_for_quantifs_scaled_log_cache.qs")


# Filter events ----
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

first_pass_regression |>
  ggplot() + theme_classic() +
  geom_point(aes(x = nb_neurons, y = nb_samples), alpha = .2) +
  xlab("Number of neurons") + ylab("Number of samples") +
  geom_hline(aes(yintercept = 120.5), color = 'darkred') +
  geom_vline(aes(xintercept = 32.5), color = 'darkred')

table(first_pass_regression$nb_samples > 120 & first_pass_regression$nb_neurons > 32)



events_to_keep <- quantifs |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  filter(nb_samples > 120,
         nb_neurons > 32) |>
  pull(event_id)


table(unique(quantifs$event_id) %in% events_to_keep)

quantifs_filtered <- quantifs |>
  filter(event_id %in% events_to_keep)


patchwork::wrap_plots(
  ggplot(first_pass_regression,
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_regression,
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ncol = 1
)



# Play with examples ----
# SE_136, SE_303 are good. 1060 particularly
(my_ev <- sample(quantifs_filtered$event_id, 1))
my_ev <- "SE_136"

quantifs_filtered |>
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


# Sparse LASSO
#scale x and remove cols that have NaNs (when there is no measurement)
x <- sf_expression |>
  select(transcript_id, sample_id, TPM) |>
  mutate(TPM = log(TPM + 1)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "TPM") |>
  column_to_rownames("sample_id") |>
  as.matrix() |>
  scale()
x <- x[,! apply(x, 2, \(col) any(is.na(col)))]


y <- quantifs_filtered[quantifs_filtered$event_id == my_ev, c("sample_id", "PSI")] |>
  column_to_rownames("sample_id") |>
  filter(!is.na(PSI)) |>
  as.matrix()

table(rownames(y) %in% rownames(x))
table(rownames(x) %in% rownames(y))

x <- x[rownames(y),]


n <- nrow(x)
train <- sample(n, round(.7*n))


cvfit <- cv.glmnet(x[train,], y[train], nfolds = 30, type.measure = "mse")
plot(cvfit)
log(cvfit$lambda.1se)
log(cvfit$lambda.min)
assess.glmnet(cvfit, newx = x[-train,], newy = y[-train], s = "lambda.min")

predict(cvfit, newx = x[-train,], s = "lambda.min") |>
  as.data.frame() |>
  as_tibble(rownames = "sample_id") |>
  add_column(measured = y[-train]) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = measured, y = lambda.min)) +
  ylab("predicted")


coef(cvfit, s = "lambda.min") |>
  as.matrix() |>
  as_tibble(rownames = "transcript_id") |>
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



#~ Randomized examples ----

my_ev <- "SE_136"
my_ev <- "SE_303"


# Sparse LASSO


x <- sf_expression |>
  select(transcript_id, sample_id, TPM) |>
  mutate(TPM = log(TPM + 1)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "TPM") |>
  column_to_rownames("sample_id") |>
  as.matrix() |>
  scale()
x <- x[,! apply(x, 2, \(col) any(is.na(col)))]
y <- quantifs[quantifs$event_id == my_ev, c("sample_id", "PSI")] |>
  column_to_rownames("sample_id") |>
  filter(!is.na(PSI)) |>
  as.matrix()

table(rownames(y) %in% rownames(x))
table(rownames(x) %in% rownames(y))

x <- x[rownames(y),]


n <- nrow(x)

# Randomize!
#all values
# x[] <- sample(x)
#by sample
x <- x[sample(nrow(x)),]

train <- sample(n, round(.7*n))



cvfit <- cv.glmnet(x[train,], y[train])
plot(cvfit)
log(cvfit$lambda.1se)
log(cvfit$lambda.min)
assess.glmnet(cvfit, newx = x[-train,], newy = y[-train], s = "lambda.min")

predict(cvfit, newx = x[-train,], s = "lambda.min") |>
  as.data.frame() |>
  as_tibble(rownames = "sample_id") |>
  add_column(measured = y[-train]) |>
  ggplot(aes(x = measured, y = lambda.min)) +
  theme_classic() +
  geom_point() +
  ylab("predicted") +
  geom_smooth(method = "lm")


coef(cvfit, s = "lambda.min") |>
  as.matrix() |>
  as_tibble(rownames = "transcript_id") |>
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

patchwork::wrap_plots(
  ggplot(first_pass_filt,
         aes(x = nb_samples, y = rsquare)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ggplot(first_pass_filt,
         aes(x = nb_samples, y = nb_coefs)) +
    theme_classic() +
    geom_point(alpha = .2) +
    geom_smooth(),
  ncol = 1
)





overlaps_sf <- first_pass_filt |>
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

overlaps_sf |>
  mutate(prop_overlapping_sf = nb_overlapping_sf/nb_known_sf) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nb_samples, y = prop_overlapping_sf), alpha = .2)
sum(overlaps_sf$nb_known_sf)
sum(overlaps_sf$nb_overlapping_sf)
126/1496
#> (no scaling) only 8% of known interactions are found in our computations
154/1750
#> 8.8% with scaling
208/1750
#> 12% on log scale

# randomize and test
randomized_overlaps <- replicate(200,{
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




# Check reproducibility ----

first_pass_regression |>
  ggplot() + theme_classic() +
  geom_point(aes(x = nb_neurons, y = nb_samples), alpha = .2) +
  xlab("Number of neurons") + ylab("Number of samples") +
  geom_hline(aes(yintercept = 178), color = 'darkred') +
  geom_vline(aes(xintercept = 45.5), color = 'darkred')

table(first_pass_regression$nb_samples > 178 &
        first_pass_regression$nb_neurons > 45)

events_to_keep2 <- quantifs |>
  filter(! is.na(PSI)) |>
  group_by(event_id) |>
  summarize(nb_samples = n(),
            nb_neurons = n_distinct(neuron_id)) |>
  filter(nb_samples > 178,
         nb_neurons > 45) |>
  pull(event_id)


table(unique(quantifs$event_id) %in% events_to_keep2)

quantifs_filtered2 <- quantifs |>
  filter(! is.na(PSI)) |>
  filter(event_id %in% events_to_keep2)





# replicated_regression <- replicate(n = 50,
#                                    expr = map(events_to_keep2, possibly(sparse_regression,
#                                                                         otherwise = list(rsquare=NA_real_,nb_coefs=NA_real_))),
#                                    simplify = FALSE)
# qs::qsave(replicated_regression, "data/intermediates/230113_replicated_regression_log_cache.qs")

replicated_regression <- qs::qread("data/intermediates/230113_replicated_regression_log_cache.qs")

replicated_rsquare <- map(replicated_regression,
                          \(replicate) {
                            map_dbl(replicate, \(rep) rep[["rsquare"]]) |>
                              set_names(events_to_keep2)
                          }) |>
  set_names(paste0("replicate_", 1:50)) |>
  as_tibble() |>
  add_column(event_id = events_to_keep2, .before = 1) |>
  pivot_longer(-event_id,
               names_to = "replicate",
               values_to = "Rsquare_adjusted")

ggplot(replicated_rsquare,
       aes(x = event_id, y = Rsquare_adjusted)) +
  theme_classic() +
  # geom_violin(fill = 'grey95') +
  geom_boxplot(fill = 'grey90') +
  geom_point() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) 


# recheck 1 value
(lm(measured ~ predicted, data = replicated_regression[[2]][[3]]$prediction_on_test) |>
    summary())$adj.r.squared

replicated_rsquare |>
  filter(replicate == "replicate_2", event_id == events_to_keep2[[3]]) |>
  pull(Rsquare_adjusted)


# also look at SF computed
replicated_regression[[1]][[1]][["coefs_sf"]]
replicated_coefs <- imap(replicated_regression,
                         \(replicate, ind) {
                           map2(replicate, events_to_keep2,
                                \(rep, .event_id) {
                                  rep[["coefs_sf"]] |>
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

sum(overlap_intersected_coefs$nb_sf_overlap)
sum(overlap_intersected_coefs$nb_sf_known)
25/482



rep_overlap <- replicate(500,
                         sum(map2_int(sample(intersected_coefs$computed), sample(intersected_coefs$known),
                                      ~ length(intersect(.x$computed_sf_gene_id, .y$sf_id))))
)

hist(rep_overlap, breaks = 30, xlab = "Number of overlapping interactions under randomization",main = NULL); abline(v = 22, col = 'red')





# Select only good fits ----
hist(replicated_rsquare$Rsquare_adjusted, breaks = 60); abline(v = .51, col = 'red')
hist(replicated_rsquare$Rsquare_adjusted, breaks = 800, xlim = c(.48,.52)); abline(v = .51, col = 'red')
table(replicated_rsquare$Rsquare_adjusted > .51)

pluck(replicated_regression[[45]][[25]],"coefs_sf", .default = tibble()) |>
  add_column(event_id = "SE111")
length(replicated_regression[[45]])
events_to_keep2[24:26]

replicated_coefs <- imap(replicated_regression,
                         \(replicate, ind) {
                           map2(replicate, events_to_keep2,
                                \(rep, .event_id) {
                                  pluck(rep, "coefs_sf", .default = tibble()) |>
                                    add_column(event_id = .event_id)
                                }) |>
                             list_rbind() |>
                             add_column(replicate = paste0("replicate_", ind))
                         }) |>
  list_rbind()





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


rep_overlap <- replicate(500,
                         sum(map2_int(sample(filt_intersected_coefs$computed), sample(filt_intersected_coefs$known),
                                      ~ length(intersect(.x$computed_sf_gene_id, .y$sf_id))))
)

hist(rep_overlap, breaks = 10, xlab = "Number of overlapping interactions under randomization",main = NULL); abline(v = 10, col = 'red')

#> not better






