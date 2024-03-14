# Prepare data for "graph_power4", so the main script consistently uses the same data



library(tidyverse)


library(wbData)


gids <- wb_load_gene_ids(289)
tx2g <- wb_load_tx2gene(289)




datadir <- "data/graph_power4/inputs/"


#~ Load data ----
quantifs_filtered <- qs::qread(file.path(datadir, "240308_preprocessed_quantifs_filtered.qs"))
sf_expression <- qs::qread(file.path(datadir, "240308_preprocessed_sf_expression.qs")) |>
  filter(transcript_id != "R07E5.14.2")


#~ Prepare data ----
message("---- Prepare data")


#~~ PSI -----
mat_psi <- quantifs_filtered |>
  select(event_id, sample_id, PSI) |>
  pivot_wider(id_cols = sample_id,
              names_from = event_id,
              values_from = PSI
  ) |>
  column_to_rownames("sample_id") |>
  as.matrix()



# filter PSI


# remove samples full of NA
mat_psi <- mat_psi[rowMeans(is.na(mat_psi)) < .4, ]


# Train/test split: note we do that BEFORE imputation
set.seed(123)
train_samples <- sample(rownames(mat_psi), size = .7*nrow(mat_psi))
test_samples <- setdiff(rownames(mat_psi), train_samples)

mat_psi_train <- mat_psi[train_samples,]




#~ SF TPM ----
mat_sf <- sf_expression |>
  mutate(logTPM = log(TPM + 1)) |>
  select(transcript_id, sample_id, logTPM) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "logTPM") |>
  column_to_rownames("sample_id") |>
  as.matrix()
mat_sf_train <- mat_sf[train_samples, ]




#~ Assemble data ----

# match rows
stopifnot(all.equal(rownames(mat_psi_train), rownames(mat_sf_train)))
mat_train <- cbind(mat_psi_train, mat_sf_train)


# finish
nb_psi <- ncol(mat_psi_train)
nb_sf <- ncol(mat_sf_train)

mat_test <- cbind(mat_sf[test_samples,], mat_psi[test_samples, ])
mat_sf_test <- mat_test[,1:nb_sf]
mat_psi_test <- mat_test[,(nb_sf+1):(nb_sf+nb_psi)]


# 5-fold cross-validation
set.seed(1)
folds <- (rep(1:5,
              each = ceiling(nrow(mat_train)/5)) |>
            sample())[1:nrow(mat_train)]


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

save_dir <- "data/graph_power4/inputs/240313_precomputed/"


qs::qsave(nb_psi, file.path(save_dir, "nb_psi.qs"))
qs::qsave(nb_sf, file.path(save_dir, "nb_sf.qs"))
qs::qsave(mat_train, file.path(save_dir, "mat_train.qs"))
qs::qsave(folds, file.path(save_dir, "folds.qs"))
qs::qsave(mat_interactions_lit, file.path(save_dir, "mat_interactions_lit.qs"))



