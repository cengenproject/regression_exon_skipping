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



library(ComplexHeatmap)

# central matrix
coefs <- as.matrix(all_coefs[rowSums(all_coefs) > 0, ])
coefs[coefs == 0] <- NA
coefs <- as.matrix(t(coefs))

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


h_full <- Heatmap(coefs[names_psi,names_x],heatmap_width = unit(0.3, "npc"), heatmap_height = unit(0.3, "npc"),
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
        left_annotation = h_psi)

pdf("data/intermediates/heatmaps/full1.pdf", width = 15, height = 15)
draw(h_full)
dev.off()


# ready-to-assemble

h_coefs <- Heatmap(coefs[names_psi,names_x],
                   rect_gp = gpar(col = "grey90", lwd = .5),
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   na_col = "white",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   col = circlize::colorRamp2(breaks = c(min(all_coefs), 0, max(all_coefs)),
                                              colors = c("red", "grey","blue")),
                   show_heatmap_legend = FALSE, left_annotation = HeatmapAnnotation(df = as.data.frame(t(psi[,names_psi]))))


h_x <- Heatmap(x[,names_x],
               show_row_names = FALSE,
               show_column_names = FALSE,
               na_col = "white",
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               col = circlize::colorRamp2(breaks = c(min(x), 0, max(x)),
                                          colors = c("magenta", "white","green4")),
               show_heatmap_legend = FALSE)



h_psi <- Heatmap(t(psi[,names_psi]),
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 na_col = "grey",
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 col = circlize::colorRamp2(breaks = c(min(psi, na.rm = TRUE), 0, max(psi, na.rm = TRUE)),
                                            colors = c("orange4", "white","purple4")),
                 show_heatmap_legend = FALSE)


fname <- "data/intermediates/heatmaps/heatmap1"
pdf(paste0(fname, "_psi.pdf"), width = 5, height = 3)
  draw(h_psi)
dev.off()
pdf(paste0(fname, "_x.pdf"), width = 3, height = 5)
  draw(h_x)
dev.off()
pdf(paste0(fname, "_coefs.pdf"), width = 3, height = 3)
  draw(h_coefs)
dev.off()






# With legend
h_coefs <- Heatmap(coefs[names_psi,names_x],
        show_row_names = FALSE,
        show_column_names = FALSE,
        na_col = "white",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = circlize::colorRamp2(breaks = c(min(all_coefs), 0, max(all_coefs)),
                                   colors = c("red", "grey","blue")),
        heatmap_legend_param = list(title = "Regression coefficient"))


h_x <- Heatmap(x[,names_x],
        show_row_names = FALSE,
        show_column_names = FALSE,
        na_col = "white",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = circlize::colorRamp2(breaks = c(min(x), 0, max(x)),
                                   colors = c("magenta", "white","green4")),
        heatmap_legend_param = list(title = "Expression Z-score"))



h_psi <- Heatmap(t(psi[,names_psi]),
        show_row_names = FALSE,
        show_column_names = FALSE,
        na_col = "grey",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = circlize::colorRamp2(breaks = c(min(psi, na.rm = TRUE), 0, max(psi, na.rm = TRUE)),
                                   colors = c("orange4", "white","purple4")),
        heatmap_legend_param = list(title = "PSI"))


 








image(t(psi[,names_psi]))

pheatmap(coefs[names_psi,names_x],
         cluster_rows = FALSE, cluster_cols = FALSE)


