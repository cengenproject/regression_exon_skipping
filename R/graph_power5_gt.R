# check results from graph_power4.R (use the cached results)
# check against ground truth connections


# Inits ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)
tx2g <- wb_load_tx2gene(281)


events_coords <- read_tsv("data/export_for_arman/221111_events_coordinates.tsv") |>
  select(event_id, gene_id)



# Get ground truth ----
all_interactions <- readr::read_tsv("../../biblio_SF/outputs/sf_targets_v3.tsv",
                                    show_col_types = FALSE)

all_interactions_by_event <- all_interactions |>
  select(SF, targets) |>
  distinct() |>
  left_join(events_coords,
            by = join_by(targets == gene_id),
            relationship = "many-to-many") |>
  filter(!is.na(event_id)) |>
  select(event_id, SF) |>
  mutate(sf_tx = wb_g2tx(SF, tx2g)) |>
  select(-SF) |>
  unnest(sf_tx) |>
  distinct() |>
  add_column(literature = TRUE)



# Functions ----
get_coefs_from_OM <- function(OM, permute = FALSE){
  
  adj <- OM[startsWith(rownames(OM), "SE_"),
            !startsWith(colnames(OM), "SE_")]
  
  if(permute){
    adj <- apply(adj, 2, sample)
    rownames(adj) <- rownames(OM[startsWith(rownames(OM), "SE_"),])
  }
  
  adj |>
    as.data.frame() |>
    rownames_to_column("event_id_percount") |>
    pivot_longer(cols = -event_id_percount,
                 names_to = "sf_id",
                 values_to = "coefficient") |>
    separate_wider_delim(event_id_percount,
                         delim = ".",
                         names = c("event_id", NA)) |>
    group_by(sf_id, event_id) |>
    summarize(coefficient = max(coefficient),
              .groups = "drop")
}




# QUIC ----
# with every rho
OM_quic_all <- qs::qread("data/intermediates/230920_cv/230921_tib_quic2.qs")



coefs_quic <- map(unique(OM_quic_all$penalty),
                     \(.penalty) {
                       
                       OM_single <- OM_quic_all |>
                         filter(penalty == .penalty) |>
                         pull(OM)
                       
                       
                       coefs <- imap(OM_single,
                                          \(.OM, .fold){
                                            map(0:5,
                                                 \(.permutation){
                                                   get_coefs_from_OM(.OM, .permutation) |>
                                                     add_column(fold = .fold,
                                                                permutation = .permutation)
                                                 }) |>
                                              list_rbind()
                                          }) |>
                         list_rbind()
                       
                       
                       coefs |>
                         group_by(permutation, sf_id, event_id) |>
                         summarize(nb_folds = sum(coefficient != 0),
                                   .groups = "drop") |>
                         left_join(all_interactions_by_event,
                                   by = c("event_id", "sf_id" = "sf_tx")) |>
                         mutate(literature = replace_na(literature, FALSE)) |>
                         add_column(penalty = .penalty,
                                    .before = 1)
                     },
                     .progress = TRUE
)

res_quic <- coefs_quic |>
  list_rbind() |>
  mutate(detected = nb_folds > 2) |>
  group_by(permutation, penalty) |>
  summarize(prop_known = 100 * sum(detected & literature)/sum(literature),
            prop_unknown = 100 * sum(detected & ! literature)/sum(!literature),
            .groups = 'drop') |>
  pivot_longer(-c(penalty, permutation),
               names_to = "literature",
               names_prefix = "prop_",
               values_to = "proportion detected")

res_quic |>
  ggplot(aes(x = penalty, y = `proportion detected`, color = permutation != 0, shape = literature)) +
  theme_classic() +
  geom_point() +
  # geom_line() +
  scale_x_log10()


res_quic |>
  group_by(penalty, literature) |>
  nest() |>
  mutate(normalized_proportion = map_dbl(data,
                                         ~ .x$`proportion detected`[.x$permutation == 0] /
                                           mean(.x$`proportion detected`[.x$permutation != 0]))) |>
  select(-data) |>
  ggplot(aes(x = penalty, y = normalized_proportion, color = literature)) +
  theme_classic() +
  geom_point() +
  geom_line() +
  scale_x_log10()








# huge ----
# with every rho
OM_huge_all <- qs::qread("data/intermediates/230920_cv/230921_tib_huge.qs")



coefs_huge_gt <- map(unique(OM_huge_all$penalty),
                     \(.penalty) {
                       
                       OM_single <- OM_huge_all |>
                         filter(penalty == .penalty) |>
                         pull(OM)
                       
                       
                       coefs <- imap(OM_single,
                                     \(.OM, .fold){
                                       map(0:5,
                                           \(.permutation){
                                             get_coefs_from_OM(.OM, .permutation) |>
                                               add_column(fold = .fold,
                                                          permutation = .permutation)
                                           }) |>
                                         list_rbind()
                                     }) |>
                         list_rbind()
                       
                       
                       coefs |>
                         group_by(permutation, sf_id, event_id) |>
                         summarize(nb_folds = sum(coefficient != 0),
                                   .groups = "drop") |>
                         left_join(all_interactions_by_event,
                                   by = c("event_id", "sf_id" = "sf_tx")) |>
                         mutate(literature = replace_na(literature, FALSE)) |>
                         add_column(penalty = .penalty,
                                    .before = 1)
                     },
                     .progress = TRUE
)

coefs_huge_gt2 <- coefs_huge_gt |>
  map2(unique(OM_huge_all$penalty),
       ~ add_column(.x, penalty = .y)) |>
  list_rbind() |>
  mutate(detected = nb_folds > 0)


coefs_huge_gt2 |>
  group_by(penalty) |>
  summarize(prop_known = 100 * sum(detected & literature)/sum(literature),
            prop_unknown = 100 * sum(detected & ! literature)/sum(!literature)) |>
  pivot_longer(-penalty,
               names_to = "literature",
               names_prefix = "prop_",
               values_to = "proportion detected") |>
  ggplot(aes(x = penalty, y = `proportion detected`, color = literature)) +
  theme_classic() +
  geom_point() +
  geom_line() +
  scale_x_log10()


coefs_huge_gt2 |>
  filter(penalty == 0.15) -> xx


table(xx$nb_folds, xx$literature)



  