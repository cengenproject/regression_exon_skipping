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







# ~~~~~~~~~~~~ ----





# QUIC ----
# with every rho
OM_quic_all <- qs::qread("data/intermediates/230920_cv/230921_tib_quic2.qs")



coefs_quic <- map(unique(OM_quic_all$penalty),
                     \(.penalty) {
                       
                       OM_single <- OM_quic_all |>
                         filter(penalty == .penalty) |>
                         pull(OM)
                       
                       #~ gather the coefficients ----
                       coefs <- imap(OM_single,
                                          \(.OM, .fold){
                                            map(0:101,
                                                 \(.permutation){
                                                   get_coefs_from_OM(.OM, .permutation) |>
                                                     add_column(fold = .fold,
                                                                permutation = .permutation)
                                                 }) |>
                                              list_rbind()
                                          }) |>
                         list_rbind()
                       
                       #~ summarize coefficients two methods ----
                       
                       coef_nonzero <- coefs |>
                         group_by(permutation, sf_id, event_id) |>
                         summarize(nb_folds_non_zero = sum(coefficient != 0),
                                   .groups = "drop") |>
                         left_join(all_interactions_by_event,
                                   by = c("event_id", "sf_id" = "sf_tx")) |>
                         mutate(literature = replace_na(literature, FALSE))
                       
                       
                       coef_signif <- coefs |>
                         group_by(sf_id, event_id, fold) |>
                         nest() |>
                         ungroup() |>
                         mutate(pval = map_dbl(data,
                                               ~ mean(.x$coefficient[.x$permutation > 0] >= 
                                                        .x$coefficient[.x$permutation == 0]))) |>
                         mutate(pval_adj = p.adjust(pval, "BH")) |>
                         select(-data) |>
                         group_by(sf_id, event_id) |>
                         summarize(nb_folds_sig = sum(pval_adj < 0.1),
                                   .groups = "drop")
                       
                       
                       #~ output proportion per penalty ----
                       prop_nonzero <- coef_nonzero |>
                         mutate(detected = nb_folds_non_zero > 2) |>
                         group_by(permutation) |>
                         summarize(prop_known = 100 * sum(detected & literature)/sum(literature),
                                   prop_unknown = 100 * sum(detected & ! literature)/sum(!literature),
                                   .groups = 'drop') |>
                         pivot_longer(-permutation,
                                      names_to = "literature",
                                      names_prefix = "prop_",
                                      values_to = "proportion_detected") |>
                         group_by(literature) |>
                         nest() |>
                         mutate(normalized_proportion_nonzero = map_dbl(data,
                                                                        ~ .x$proportion_detected[.x$permutation == 0] /
                                                                          mean(.x$proportion_detected[.x$permutation != 0]))) |>
                         select(-data)
                       
                       prop_signif <- coef_signif |>
                         left_join(all_interactions_by_event,
                                   by = c("event_id", "sf_id" = "sf_tx")) |>
                         mutate(literature = replace_na(literature, FALSE),
                                literature = if_else(literature, "known", "unknown")) |>
                         group_by(literature) |>
                         summarize(proportion_detected = 100 * mean(nb_folds_sig > 2))
                       
                       #~ gather and export ----
                       left_join(prop_nonzero, prop_signif,
                                              by = "literature") |>
                         add_column(penalty = .penalty,
                                    .before = 1)
                       
                     },
                     .progress = TRUE
) |>
  list_rbind()

# qs::qsave(coefs_quic, "data/intermediates/231002_permutations/231003_coefs_quic.qs")




#~ plot ----


coefs_quic |>
  ggplot(aes(x = penalty, y = normalized_proportion_nonzero, color = literature)) +
  theme_classic() +
  geom_point() +
  geom_line() +
  scale_x_log10()


#~ ground truth ----
OM <- (OM_quic_all |>
  filter(penalty == .1,
         fold == 1) |>
  pull(OM))[[1]]
coefs <- get_coefs_from_OM(OM)

dim(coefs)
coefs[1:3,]

# huge ----
# with every rho
OM_huge_all <- qs::qread("data/intermediates/230920_cv/230921_tib_huge.qs")



coefs_huge <- map(unique(OM_huge_all$penalty),
                     \(.penalty) {
                       
                       OM_single <- OM_huge_all |>
                         filter(penalty == .penalty) |>
                         pull(OM)
                       
                       
                       coefs <- imap(OM_single,
                                     \(.OM, .fold){
                                       map(0:51,
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
# qs::qsave(coefs_huge, "data/intermediates/231002_permutations/231002_coefs_huge.qs")


res_huge <- coefs_huge |>
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

res_huge |>
  ggplot(aes(x = penalty, y = `proportion detected`, color = permutation != 0, shape = literature)) +
  theme_classic() +
  geom_point() +
  # geom_line() +
  scale_x_log10()


res_huge |>
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









# ARACNE ----
# with every rho
OM_arac_all <- qs::qread("data/intermediates/230920_cv/230922_tib_arac_raw.qs")
          
#~ gather the coefficients ----
OM_single <- OM_arac_all |>
  pull(OM)

coef_arac <- imap(OM_single,
              \(.OM, .fold){
                map(0:51,
                    \(.permutation){
                      get_coefs_from_OM(.OM, .permutation) |>
                        add_column(fold = .fold,
                                   permutation = .permutation)
                    }) |>
                  list_rbind()
              },
              .progress = TRUE) |>
  list_rbind()
# qs::qsave(coef_arac, "data/intermediates/231002_permutations/231003_coef_arac.qs")


#~ summarize coefficients two methods ----
coef_signif <- coef_arac |>
  group_by(sf_id, event_id, fold) |>
  nest() |>
  ungroup() |>
  mutate(pval = map_dbl(data,
                        ~ mean(.x$coefficient[.x$permutation > 0] >= 
                                 .x$coefficient[.x$permutation == 0]))) |>
  mutate(pval_adj = p.adjust(pval, "BH")) |>
  select(-data) |>
  group_by(sf_id, event_id) |>
  summarize(nb_folds_sig = sum(pval_adj < 0.1),
            .groups = "drop")

coef_nonzero <- coef_arac |>
  group_by(permutation, sf_id, event_id) |>
  summarize(nb_folds_non_zero = sum(coefficient != 0),
            .groups = "drop") |>
  left_join(all_interactions_by_event,
            by = c("event_id", "sf_id" = "sf_tx")) |>
  mutate(literature = replace_na(literature, FALSE))


#~ output proportion per penalty ----
prop_nonzero <- coef_nonzero |>
  mutate(detected = nb_folds_non_zero > 2) |>
  group_by(permutation) |>
  summarize(prop_known = 100 * sum(detected & literature)/sum(literature),
            prop_unknown = 100 * sum(detected & ! literature)/sum(!literature),
            .groups = 'drop') |>
  pivot_longer(-permutation,
               names_to = "literature",
               names_prefix = "prop_",
               values_to = "proportion_detected") |>
  group_by(literature) |>
  nest() |>
  mutate(normalized_proportion_nonzero = map_dbl(data,
                                                 ~ .x$proportion_detected[.x$permutation == 0] /
                                                   mean(.x$proportion_detected[.x$permutation != 0]))) |>
  select(-data)

prop_signif <- coef_signif |>
  left_join(all_interactions_by_event,
            by = c("event_id", "sf_id" = "sf_tx")) |>
  mutate(literature = replace_na(literature, FALSE),
         literature = if_else(literature, "known", "unknown")) |>
  group_by(literature) |>
  summarize(proportion_detected = 100 * mean(nb_folds_sig > 2))

prop_arac <- left_join(prop_nonzero, prop_signif,
                       by = "literature")


#~ plot ----

prop_arac |>
  pivot_longer(-literature,
               names_to = "metric",
               values_to = "proportion") |>
  ggplot() +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_point(aes(x = literature, y = proportion))












# DPM ----
# with every rho
OM_dpm_all <- qs::qread("data/intermediates/230920_cv/230922_tib_dpm.qs")



coefs_dpm <- map(unique(OM_dpm_all$penalty),
                  \(.penalty) {
                    
                    OM_single <- OM_dpm_all |>
                      filter(penalty == .penalty) |>
                      pull(OM)
                    
                    
                    coefs <- imap(OM_single,
                                  \(.OM, .fold){
                                    map(0:51,
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
# qs::qsave(coefs_dpm, "data/intermediates/231002_permutations/231002_coefs_dpm.qs")


res_dpm <- coefs_dpm |>
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

res_dpm |>
  ggplot(aes(x = penalty, y = `proportion detected`, color = permutation != 0, shape = literature)) +
  theme_classic() +
  geom_point() +
  # geom_line() +
  scale_x_log10()


res_huge |>
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





  