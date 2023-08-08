# Reproduce the scQuint paper fig 11
# https://elifesciences.org/articles/73520
# 
# Methods:
# We filtered target exon skipping events to those defined in at least 95% of the replicates,
# and those having a PSI standard deviation of at least 0.2. We used log-transformed 
# normalized expression and PSI of alternative splicing events as input features.
# We chose to keep the PSI of only one intron per intron group to avoid the presence of highly
# correlated features and improve clarity, even if some information from non-binary events is lost. 
# Input features were filtered to those having standard deviation of at least 0.05, and then standardized. 
# A lasso Dirichlet-Multinomial GLM was fit to the data (in this instance, the model reduces to a
# Beta-Binomial because skipped exons are binary events), with the sparsity penalty selected via 
# cross-validation. As a first approach, we fit a regular lasso linear regression model on PSI 
# instead of raw counts, resulting in roughly similar patterns in the coefficients. Figure 11c
# shows the coefficients of the lasso Dirichlet-Multinomial model for the top 30 targets with the
# highest variance explained by the regular lasso model, all above 68%.



library(tidyverse)
library(glmnet)
library(wbData)

gids <- wb_load_gene_ids(281)
tx2g <- wb_load_tx2gene(281)

quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs")

events_coords <- read_tsv("data/export_for_arman/221111_events_coordinates.tsv") |>
  select(event_id, gene_id)


list_sfs <- unique(sf_expression$gene_id)


# Input features ----
# normalized (i.e. TPM), log, expression; and PSI of events (from SFs)
x_sf <- sf_expression |>
  mutate(TPM = log1p(TPM)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "TPM") |>
  column_to_rownames("sample_id") |>
  as.matrix()
x_psi <- quantifs_filtered |>
  left_join(events_coords, by = "event_id") |>
  filter(gene_id %in% list_sfs) |>
  mutate(transcript_id = paste(gene_id, event_id, sep = "_")) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "PSI") |>
  column_to_rownames("sample_id") |>
  as.matrix()

x_psi <- x_psi[,colMeans(is.na(x_psi)) < 0.2]


# filter
sds_sf <- matrixStats::colSds(x_sf)
sds_psi <- matrixStats::colSds(x_psi, na.rm = TRUE)

hist(sds_sf, breaks = 50); abline(v = 0.7, col ='darkred', lty = 'dashed')
hist(sds_psi, breaks = 50)

x_psi <- impute::impute.knn(x_psi)$data

x <- cbind(x_sf[,sds_sf > 0.7], x_psi[, sds_psi > 0.1]) |>
  scale()


dim(x)
image(x)






# targets

keep_events <- quantifs_filtered |>
  complete(event_id, sample_id) |>
  group_by(event_id) |>
  summarize(prop_na = mean(is.na(PSI)),
            psi_sd = sd(PSI, na.rm = TRUE)) |>
  filter(prop_na < .4,
         psi_sd > .2) |>
  pull(event_id)

# binarize PSI (to apply binomial)
hist(quantifs_filtered$PSI, breaks = 50)

y <- quantifs_filtered |>
  filter(event_id %in% keep_events) |>
  mutate(PSI = 1L*(PSI > 0.5)) |>
  pivot_wider(id_cols = sample_id,
              names_from = "event_id",
              values_from = PSI) |>
  column_to_rownames("sample_id")



# filter out samples full of NA
keep_samples <- rownames(y)[rowMeans(is.na(y)) < .2]
x <- x[keep_samples, ]
y <- y[keep_samples, ]


# Fit!

mods <- map(y,
            \(col){
              glmnet(x = x[!is.na(col),], y = col[!is.na(col)], family = "binomial", intercept = FALSE)
            })

all_coefs <- do.call(cbind, map(mods, coef, s = 0.05))
colnames(all_coefs) <- names(mods)

# coefs <- all_coefs[rowSums(all_coefs) > 0, ]
# image(coefs, xlab = NULL, ylab = NULL, main = NULL)
# 
# image(x)



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
  filter(sf_tx %in% colnames(x))


reported_inter <- all_interactions_by_event |>
  mutate(reported = 1) |>
  bind_rows(tibble(event_id = setdiff(colnames(all_coefs), all_interactions_by_event$event_id),
                   sf_tx = NA,
                   reported = 0)) |>
  bind_rows(
    tibble(event_id = NA,
           sf_tx = setdiff(colnames(x), all_interactions_by_event$sf_tx),
           reported = 0)) |>
  pivot_wider(id_cols = event_id,
              names_from = "sf_tx",
              values_from = "reported") |>
  filter(!is.na(event_id)) |>
  column_to_rownames("event_id") |>
  as.matrix() |>
  (\(mat) {mat[is.na(mat)] <- 0; mat})()

table(reported_inter)
dim(reported_inter)
reported_inter[1:3,1:3]




# Accuracy ----
coefs <- t(as.matrix(all_coefs[rowSums(all_coefs) > 0, ]))
TPR <- sum((coefs !=0) * reported_inter, na.rm = TRUE)/sum(reported_inter)





library(ComplexHeatmap)

# central matrix
coefs[coefs == 0] <- NA

# left matrix
psi <- quantifs_filtered |>
  filter(event_id %in% keep_events,
         sample_id %in% keep_samples) |>
  pivot_wider(id_cols = sample_id,
              names_from = "event_id",
              values_from = PSI) |>
  column_to_rownames("sample_id") |>
  as.matrix()




# prepare order
hc_x <- hclust(dist(t(x)))
hc_psi <- hclust(dist(t(psi)))

names_x <- hc_x$labels[hc_x$order] |> intersect(colnames(coefs))
names_psi <- hc_psi$labels[hc_psi$order] |> intersect(rownames(coefs))

stopifnot(all(names_psi %in% rownames(reported_inter)))
stopifnot(all(names_x %in% colnames(reported_inter)))
reported_inter <- reported_inter[names_psi, names_x]
reported_inter_bool <- reported_inter == 1
# dim(reported_inter_bool)
# table(reported_inter_bool)
# reported_inter_bool[1:3,1:3]

h_x <- HeatmapAnnotation(x = t(x[,names_x]),
                         which = "column",
                         show_annotation_name = FALSE,
                         na_col = "white",
                         simple_anno_size = unit(1, "mm"),
                         col = list(x = circlize::colorRamp2(breaks = c(min(x), 0, max(x)),
                                                             colors = c("magenta", "white","green4"))))

h_psi <- HeatmapAnnotation(psi = t(psi[,names_psi]),
                           which = "row",
                           show_annotation_name = FALSE,
                           na_col = "grey",
                           simple_anno_size = unit(1, "mm"),
                           col = list(psi = circlize::colorRamp2(breaks = c(min(psi, na.rm = TRUE), 0, max(psi, na.rm = TRUE)),
                                                                 colors = c("orange4", "white","purple4"))))

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(if_else(reported_inter_bool[i,j],"X",""), x, y, gp = gpar(fontsize = 8, family = "bold"))
}


h_full <- Heatmap(coefs[names_psi,names_x],
                  heatmap_width = unit(0.3, "npc"),
                  heatmap_height = unit(0.3, "npc"),
                  rect_gp = gpar(col = "grey90", lwd = .5),
                  column_names_side = "top",
                  row_names_side = "right",
                  na_col = "white",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  col = circlize::colorRamp2(breaks = c(min(all_coefs), 0, max(all_coefs)),
                                             colors = c("red", "grey","blue")),
                  # show_heatmap_legend = FALSE,
                  bottom_annotation = h_x,
                  left_annotation = h_psi,
                  cell_fun = cell_fun)

draw(h_full)

pdf("data/intermediates/heatmaps/s0-05.pdf", width = 14.4, height = 9.6)
draw(h_full)
dev.off()





pdf("data/intermediates/heatmaps/full2.pdf", width = 15, height = 15)
draw(h_full)
dev.off()





# Accuracy: range of sparsity ----

get_tpr <- function(s){
  all_coefs <- do.call(cbind, map(mods, coef, s = s))
  colnames(all_coefs) <- names(mods)
  
  
  reported_inter <- all_interactions_by_event |>
    mutate(reported = 1) |>
    bind_rows(tibble(event_id = setdiff(colnames(all_coefs), all_interactions_by_event$event_id),
                     sf_tx = NA,
                     reported = 0)) |>
    bind_rows(
      tibble(event_id = NA,
             sf_tx = setdiff(colnames(x), all_interactions_by_event$sf_tx),
             reported = 0)) |>
    pivot_wider(id_cols = event_id,
                names_from = "sf_tx",
                values_from = "reported") |>
    filter(!is.na(event_id)) |>
    column_to_rownames("event_id") |>
    as.matrix() |>
    (\(mat) {mat[is.na(mat)] <- 0; mat})()
  
  
  coefs <- t(as.matrix(all_coefs[rowSums(all_coefs) > 0, ]))
  
  reported_inter <- reported_inter[rownames(coefs), colnames(coefs)]
  
  if(sum(reported_inter) == 0) return(0)
  
  TPR <- sum((coefs !=0) * reported_inter, na.rm = TRUE)/sum(reported_inter)
  # list(sum((coefs !=0) * reported_inter, na.rm = TRUE),sum(reported_inter))
  TPR
}

get_tpr(0.05)

tpr_range <- tibble(s = mods[[1]]$lambda,
                    TPR = map_dbl(s, get_tpr))


ggplot(tpr_range) +
  theme_classic() +
  geom_point(aes(x = s, y = TPR))










