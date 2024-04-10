# Prepare data for "graph_power4", so the main script consistently uses the same data



library(tidyverse)


library(wbData)


gids <- wb_load_gene_ids(289)
tx2g <- wb_load_tx2gene(289)




datadir <- "data/graph_power4/inputs/"


#~ Load data ----
quantifs_filtered <- qs::qread(file.path(datadir, "240410_preprocessed_quantifs_filtered.qs"))
sf_expression <- qs::qread(file.path(datadir, "240410_preprocessed_sf_expression.qs")) |>
  filter(transcript_id != "R07E5.14.2")


#~ Prepare data ----
message("---- Prepare data")


#~~ PSI -----
mat_se_psi <- quantifs_filtered |>
  select(event_id, sample_id, PSI) |>
  pivot_wider(id_cols = sample_id,
              names_from = event_id,
              values_from = PSI
  ) |>
  column_to_rownames("sample_id") |>
  as.matrix()

mat_se_cnt <- quantifs_filtered |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads)) |>
  select(event_id, sample_id, Nincl, Nexcl) |>
  pivot_wider(id_cols = sample_id,
              names_from = event_id,
              values_from = c(Nincl, Nexcl),
              names_vary = "slowest",
              names_glue = "{event_id}.{.value}"
  ) |>
  column_to_rownames("sample_id") |>
  as.matrix()



# filter PSI


# remove samples full of NA
prop_missing_per_sample <- rowMeans(is.na(mat_se_psi))
mat_se_psi <- mat_se_psi[prop_missing_per_sample < .4, ]
mat_se_cnt <- mat_se_cnt[prop_missing_per_sample < .4, ]


#~ SF TPM ----
mat_sf <- sf_expression |>
  mutate(logTPM = log(TPM + 1)) |>
  select(transcript_id, sample_id, logTPM) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "logTPM") |>
  column_to_rownames("sample_id") |>
  as.matrix()



# Train/test split
set.seed(123)
all_samples <- rownames(mat_se_psi)
train_samples <- sample(all_samples, size = .7*length(all_samples))
test_samples <- setdiff(all_samples, train_samples)

nb_se <- ncol(mat_se_psi)
nb_sf <- ncol(mat_sf)









#~ Assemble data ----


mat_train_psi <- cbind(mat_se_psi[train_samples,], mat_sf[train_samples, ])
mat_train_cnt <- cbind(mat_se_cnt[train_samples,], mat_sf[train_samples, ])

mat_test_psi <- cbind(mat_se_psi[test_samples,], mat_sf[test_samples, ])
mat_test_cnt <- cbind(mat_se_cnt[test_samples,], mat_sf[test_samples, ])




# 5-fold cross-validation
set.seed(123)
folds <- (rep(1:5,
              each = ceiling(nrow(mat_train_psi)/5)) |>
            sample())[seq_along(train_samples)]



# Get Ground truth ----


events_coords <- read_tsv(file.path(datadir, "240308_events_coordinates.tsv"),
                          show_col_types = FALSE) |>
  select(event_id, gene_id)

all_interactions <- read_tsv(file.path(datadir, "sf_targets_v4.tsv"),
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

mat_interactions_lit <- all_interactions_by_event |>
  pivot_wider(id_cols = event_id,
              names_from = sf_tx,
              values_from = literature) |>
  column_to_rownames("event_id") |>
  as.matrix()

mat_interactions_lit[is.na(mat_interactions_lit)] <- FALSE

# # checks
# table(colnames(mat_interactions_lit) %in% colnames(adj_mat))
# table(colnames(adj_mat) %in% colnames(mat_interactions_lit))
# 
# table(rownames(mat_interactions_lit) %in% rownames(adj_mat))
# table(rownames(adj_mat) %in% rownames(mat_interactions_lit))
# 
# not_found <- rownames(adj_mat)[!rownames(adj_mat) %in% rownames(mat_interactions_lit)]
# 
# events_coords |> filter(!event_id %in% not_found) |>
#   mutate(gene_name = i2s(gene_id, gids))



# Save preprocessed data ----

save_dir <- "data/graph_power4/inputs/240410_precomputed/"


qs::qsave(nb_se, file.path(save_dir, "nb_se.qs"))
qs::qsave(nb_sf, file.path(save_dir, "nb_sf.qs"))
qs::qsave(mat_train_psi, file.path(save_dir, "mat_train_psi.qs"))
qs::qsave(mat_train_cnt, file.path(save_dir, "mat_train_cnt.qs"))
qs::qsave(folds, file.path(save_dir, "folds.qs"))
qs::qsave(mat_interactions_lit, file.path(save_dir, "mat_interactions_lit.qs"))


qs::qsave(mat_test_psi, file.path(save_dir, "mat_test_psi.qs"))
qs::qsave(mat_test_cnt, file.path(save_dir, "mat_test_cnt.qs"))


