# Illustrations

library(tidyverse)

params <- list(
  date = "240314g",
  exonsInput = "PSI",
  transformation = "npnshrink",
  imputation = "knn",
  permutations = "0",
  penalties = "c(10, 2, 1, .5)",
  # penalties = "c(10, 5, 2, 1, .7, .5, .4, .3, .2, .1)",,
  knn_k = 10,
  algo = "glasso"
)

source("R/loss_functions.R")
source("R/functions_steps.R")

datadir <- "data/graph_power4/inputs/240314_precomputed/"

nb_se <- qs::qread(file.path(datadir, "nb_se.qs"))
nb_sf <- qs::qread(file.path(datadir, "nb_sf.qs"))
folds <- qs::qread(file.path(datadir, "folds.qs"))
mat_interactions_lit <- qs::qread(file.path(datadir, "mat_interactions_lit.qs"))


mat_train <- qs::qread(file.path(datadir, "mat_train_psi.qs"))



nb_se <- switch (params$exonsInput,
                 PSI = nb_se,
                 counts = 2*nb_se
)

mat_train |> dim()

pheatmap::pheatmap(mat_train,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   scale = "column",
                   annotation_col = data.frame(type = rep(c("SE", "SF"),
                                                          times = c(nb_se, nb_sf)),
                                               row.names = colnames(mat_train)))


# Data matrix ----
pheatmap::pheatmap(mat_train[,1:nb_se]-0.5,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   color = colorRampPalette(c("orange3", "brown3"))(10),
                   show_rownames = FALSE, show_colnames = FALSE,
                   legend = FALSE,
                   filename = "presentations/240403_fig/mat_se.png",
                   width = .5,
                   height = 1)

pheatmap::pheatmap(log1p(mat_train[,(nb_se+1):(nb_se+nb_sf)]),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   color = colorRampPalette(c("cyan4", "green4"))(10),
                   show_rownames = FALSE, show_colnames = FALSE,
                   legend = FALSE,
                   filename = "presentations/240403_fig/mat_sf.png",
                   width = 2,
                   height = 1)

# Cov matrix ----

# run the beginning of graph_power4
cov_mat <- res_quic1$S_train_t[[1]]

pheatmap::pheatmap(cov_mat,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE)



