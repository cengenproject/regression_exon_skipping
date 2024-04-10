# Network analysis
#
# From cluster

library(tidyverse)
library(igraph)

library(wbData)


tx2gt <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)

main <- qs::qread("data/graph_power4/from_cluster/archive/final_save/240404_glasso_PSI_npnshrink_knn_k10_2_16.qs")|>
  filter(penalty == 0.3, permutation == 0)


adj <- main$adj[[1]]




#~ Plot network ----
#remove unconnected
adj2 <- adj[,colSums(abs(adj)) > 0]
adj2 <- adj2[rowSums(abs(adj2)) > 0,]

nb_sf <- ncol(adj2)
nb_se <- nrow(adj2)

adjm <- rbind(
  cbind(matrix(0, nrow = nb_sf, ncol = nb_sf), t(adj2)),
  cbind(adj2, matrix(0, nrow = nb_se, ncol = nb_se))
)
adjm <- adjm != 0
dimnames(adjm) <- list(c(colnames(adj2), rownames(adj2)),
                       c(colnames(adj2), rownames(adj2)))

gr <- graph_from_adjacency_matrix(adjm,
                                  mode = "directed")



V(gr)$type <- rep(c("SF", "SE"), times = c(nb_sf, nb_se))
V(gr)$color <- rep(c("red4", "blue4"), times = c(nb_sf, nb_se))

gr

plot(gr,
     vertex.label = NA,
     vertex.size = 4,
     edge.width = 0.1,
     edge.arrow.mode = 0,
     layout = layout.graphopt(gr, spring.length = 100, spring.constant = .01))


plot(gr,
     vertex.label = NA,
     vertex.size = 4,
     edge.width = 0.1,
     edge.arrow.mode = 0,
     layout = layout.star)


gr

degree(gr) |> table()


#~ Network ----



dim(adj)
adj[1:3,1:3]
image(adj)

pheatmap::pheatmap(adj)



hist(colSums(adj))



connectivity <- colSums(adj != 0)

tab <- table(connectivity)
if(length(tab) == 1L) return(NA_real_)

k <- as.numeric(names(tab))[-1]
p_k <- as.numeric(tab)[-1]

plot((k), (p_k))

cor(log10(k), log10(p_k))^2




## special cases



tx_unc75 <- wb_g2tx(s2i("unc-75", gids), tx2gt, simplify = TRUE)
adj[, tx_unc75] |> colSums()


colSums(abs(adj)) |> hist()
abs(adj)[, tx_unc75] |> colSums()

abline(v = abs(adj)[, tx_unc75] |> colSums())

names(colSums(abs(adj)))[colSums(abs(adj)) > .5] |>
  wb_tx2g(tx2gt) |> i2s(gids)


