# Network analysis
#
# From cluster

# Inits ----

library(tidyverse)
library(igraph)

library(ggraph)
library(tidygraph)

library(wbData)


tx2gt <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)


se_coords <- read_tsv("data/graph_power4/inputs/240308_events_coordinates.tsv")

get_target_id <- function(event_id){
  res <- se_coords$gene_id[match(event_id, se_coords$event_id, incomparables = NA)]
  if (any(is.na(res))) {
    warning("get_target_id: ", sum(is.na(res)), " event IDs could not be converted. NA are returned.")
  }
  res
}

# main <- qs::qread("data/graph_power4/from_cluster/240430_final/240426_final_glasso_PSI_npntrunc_knn_k4_2_4.qs")|>
#   filter(penalty == 0.25, permutation == 0)

main <- qs::qread("data/graph_power4/from_cluster/240426_final/240426_final_glasso_PSI_npntrunc_knn_k4_2_16.qs")|>
  filter(penalty == 0.25, permutation == 0)

adj <- main$adj[[1]]





# degrees by gene ----

# count only once per gene - SE
degree_by_gene <- as.data.frame(adj) |>
  rownames_to_column("SE") |>
  pivot_longer(-SE,
               names_to = "tx_id",
               values_to = "degree") |>
  mutate(gene_id = wb_tx2g(tx_id, tx2gt, warn_missing = TRUE),
         gene_name = i2s(gene_id, gids, warn_missing = TRUE)) |>
  mutate(degree = degree != 0) |>
  select(- tx_id) |>
  distinct() |>
  summarize(degree = sum(degree),
            n = n(),
            .by = c(gene_id, gene_name)) |>
  arrange(desc(degree))

ggplot(degree_by_gene) +
  theme_classic() +
  geom_histogram(aes(x = degree), color = 'white', bins = 85) +
  xlab("Out degree") + ylab("Number of SF")



# Check which associations ----

adj_tbl <- as.data.frame(adj) |>
  rownames_to_column("event_id") |>
  pivot_longer(-event_id,
               names_to = "sf_tx",
               values_to = "degree") |>
  mutate(sf_id = wb_tx2g(sf_tx, tx2gt, warn_missing = TRUE),
         sf_name = i2s(sf_id, gids, warn_missing = TRUE),
         target_id = get_target_id(event_id),
         target_name = i2s(target_id, gids, warn_missing = TRUE))

biggest <- function(x){
  x[which.max(abs(x))]
}

adj_tbl_gene <- adj_tbl |>
  summarize(degree = biggest(degree),
            .by = c(sf_id, sf_name,target_id,target_name))

# ground truth
ground_truth <- tibble(sf_name = c("unc-75","exc-7","prp-40","mbl-1","mec-8","msi-1",
                                     "prp-40","unc-75","unc-75","unc-75","unc-75","unc-75",
                                     "ptb-1","asd-1","hrpa-1","hrpf-1","exc-7","rsp-8"),
                         target_name = c("unc-16","unc-16","unc-16","sad-1","sad-1",
                                         "sad-1","sad-1","ret-1","C07A12.7","cdgs-1",
                                         "unc-16","daf-2","daf-2","daf-2","daf-2","daf-2",
                                         "daf-2","daf-2"),
                         is_ground_truth = TRUE)

ground_truth |>
  left_join(adj_tbl_gene,
            by = c("sf_name", "target_name"))

  

adj_tbl_gene |>
  left_join(ground_truth,
            by = c("sf_name", "target_name")) |>
  replace_na(list(is_ground_truth = FALSE)) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = is_ground_truth, y = abs(degree)), alpha = .2)


adj_tbl_gene |>
  arrange(desc(abs(degree)))

hist(adj_tbl_gene$degree[adj_tbl_gene$degree != 0], breaks = 150)

#~ examples ----
adj_tbl_gene |>
  filter(sf_name == "unc-75") |>
  pull(degree) |> hist(breaks = 50)

adj_tbl_gene |>
  filter(sf_name == "exc-7") |>
  pull(degree) |> hist(breaks = 50)

adj_tbl_gene |>
  filter(sf_name == "msi-1") |>
  pull(degree) |> hist(breaks = 50)

adj_tbl_gene |>
  filter(sf_name == "asd-1") |> #filter(degree != 0)
  pull(degree) |> hist(breaks = 50)

adj_tbl_gene |>
  filter(abs(degree) > 0.05) |>
  summarize(prop_positive = mean(degree > 0),
            n = n(),
            .by = "sf_name") |>
  arrange(desc(n)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = n, y = prop_positive))


#~~ by target ----
adj_tbl_gene |>
  filter(target_name == "unc-16") |>
  filter(degree != 0)

adj_tbl_gene |>
  filter(target_name == "lin-10") |>
  filter(degree != 0)

adj_tbl_gene |>
  filter(target_name == "ret-1") |>
  filter(degree != 0)


adj_tbl_gene |>
  filter(target_name == "C07A12.7") |>
  filter(degree != 0)



adj_tbl_gene |>
  filter(target_name == "cdgs-1") |>
  filter(degree != 0)



adj_tbl_gene |>
  filter(target_name == "daf-2") |>
  filter(degree != 0)
adj_tbl_gene |>
  filter(target_name == "daf-2") |>
  filter(sf_name %in% c("unc-75", "ptb-1", "asd-1", "hrpa-1", "hrpf-1", "exc-7", "rsp-8"))



adj_tbl_gene |>
  filter(target_name == "unc-13") |>
  filter(degree != 0)
adj_tbl_gene |>
  filter(target_name == "unc-13") |>
  filter(sf_name %in% c("exc-7", "mbl-1"))

adj_tbl_gene |>
  filter(target_name == "sad-1") |>
  filter(degree != 0)
adj_tbl_gene |>
  filter(target_name == "sad-1") |>
  filter(sf_name %in% c("mec-8", "mbl-1", "msi-1", "prp-40"))








# microexons ----

se_coords
suppa_microexons <- read_csv("../suppa_events/data/outs/240426_microexons.csv")
microexons <- suppa_microexons |>
  separate_wider_regex(event_id,
                       patterns = c(
                         gene_id = "^WBGene[0-9]{8}", ";",
                         "SE:",
                         chr = "[IVX]{1,3}", ":",
                         c1 = "[0-9]+", "-",
                         c2 = "[0-9]+", ":",
                         c3 = "[0-9]+", "-",
                         c4 = "[0-9]+", ":",
                         strand = "[+-]$"
                       )) |>
  mutate(across(c1:c4, as.integer)) |>
  mutate(
    # Upstream intron
    upstream_intron_start = c1 + 1,
    upstream_intron_end = c2 - 1,
    upstream_intron_length = upstream_intron_end - upstream_intron_start + 1,
    
    # Downstream intron
    downstream_intron_start = c3 + 1,
    downstream_intron_end = c4 - 1,
    downstream_intron_length = downstream_intron_end - downstream_intron_start + 1,
    
    # Exon
    exon_start = c2,
    exon_end = c3,
    exon_length = exon_end - exon_start + 1
  ) |>
  select(-(c1:c4))



gene_coords <- wb_load_gene_coords(289) |>
  filter(gene_id %in% microexons$gene_id) |>
  select(gene_id, gene_name = name, gene_start = start, gene_end = end, strand)


prep_se_for_microex <- se_coords |>
  filter(startsWith(event_id, "SE"),
         gene_id %in% microexons$gene_id) |>
  left_join(gene_coords,
            by = "gene_id") |>
  mutate(
    upstream_intron_start = if_else(strand == "-",
                                    gene_end - intron_end + 1,
                                    gene_start + intron_start-1),
    
    upstream_intron_end =if_else(strand == "-",
                                 gene_end - exon_end ,
                                 gene_start + exon_start-2),
    
    downstream_intron_start = if_else(strand == "-",
                                      gene_end - exon_start +2,
                                      gene_start + exon_end),
    
    downstream_intron_end = if_else(strand == "-",
                                    gene_end - intron_start + 1,
                                    gene_start + intron_end-1)
  ) |>
  select(event_id, upstream_intron_start, upstream_intron_end, downstream_intron_start, downstream_intron_end)


se_in_microex <- microexons |>
  left_join(
    prep_se_for_microex,
    by = c("upstream_intron_start", "upstream_intron_end", "downstream_intron_start",
           "downstream_intron_end")
  ) |>
  pull(event_id)


size_tot <- adj_tbl |>
  summarize(degree = biggest(degree),
            .by = c(sf_id, sf_name,target_id,target_name, event_id)) |>
  filter(degree != 0) |>
  mutate(microexon = (event_id %in% se_in_microex)) |>
  nrow()

size_micro <- adj_tbl |>
  summarize(degree = biggest(degree),
            .by = c(sf_id, sf_name,target_id,target_name, event_id)) |>
  filter(degree != 0) |>
  mutate(microexon = (event_id %in% se_in_microex)) |>
  filter(microexon) |>
  nrow()




test_res <- adj_tbl |>
  summarize(degree = biggest(degree),
            .by = c(sf_id, sf_name,target_id,target_name, event_id)) |>
  filter(degree != 0) |>
  mutate(microexon = (event_id %in% se_in_microex)) |>
  summarize(n_tot = n(),
            n_microexon = sum(microexon),
            .by = c("sf_id", "sf_name")) |>
  mutate(p = phyper(n_microexon, n_tot,
                     size_tot - n_tot,size_micro, lower.tail = FALSE),
         padj = p.adjust(p, "fdr")) |>
  arrange(padj)

test_res |>
  filter(padj < .1)
test_res |>
  arrange(desc(n_microexon))


#~ SF GO ----

# adj_tbl_gene$sf_id |>
#   unique() |>
  # writeLines("presentations/240308_figures/background.txt")

adj_tbl_gene |>
  filter(abs(degree) > 0.05) |>
  summarize(prop_positive = mean(degree > 0),
            n = n(),
            .by = "sf_id") |>
  filter(n > 4) |>
  pull(sf_id) |>
  clipr::write_clip()

adj_tbl_gene |>
  filter(degree != 0) |>
  summarize(n = n(),
            .by = "sf_id") |>
  filter(n > 4) |>
  pull(sf_id) |>
  clipr::write_clip()

adj_tbl_gene |>
  filter(degree != 0) |>
  summarize(prop_positive = mean(degree > 0),
            n = n(),
            .by = "sf_name") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = n, y = prop_positive))


# no GO enrichment of SF identity: maybe because background list is already all SF




#~ Compare "bad" ground truth ----

all_interactions <- read_tsv(file.path("data/graph_power4/inputs/", "sf_targets_v4.tsv"),
                             show_col_types = FALSE)

int_gt <- all_interactions |>
  select(SF, targets) |>
  distinct() |>
  rename(sf_id = SF,
         target_id = targets) |>
  mutate(is_ground_truth = TRUE)

int_gt |>
  left_join(adj_tbl_gene,
            by = c("sf_id", "target_id"))



adj_tbl_gene |>
  left_join(int_gt,
            by = c("sf_id", "target_id")) |>
  replace_na(list(is_ground_truth = FALSE)) |>
  filter(degree != 0) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = is_ground_truth, y = abs(degree)), alpha = .2)


adj_tbl_gene |>
  summarize(nonzero = sum(degree > 0),
            n = n(),
            .by = c(sf_id, sf_name)) |>
  arrange(desc(nonzero))



  
# ...


adj_gene <- adj_abs

sf_centrality <- colSums(adj != 0)

sort(sf_centrality, decreasing = TRUE) |> head()

enframe(sf_centrality,
        name = "tx_id",
        value = "degree") |>
  mutate(gene_id = wb_tx2g(tx_id, tx2gt, warn_missing = TRUE),
         gene_name = i2s(gene_id, gids, warn_missing = TRUE)) |>
  summarize(degree = sum(degree),
            .by = c(gene_id, gene_name)) |>
  arrange(desc(degree))

hist(sf_centrality, breaks = 150)

names(sf_centrality)[sf_centrality > 10] |>
  wb_tx2g(tx2gt) |> i2s(gids)




# plot ----
tbl <- as.data.frame(adj2) |>
  rownames_to_column("SE") |>
  pivot_longer(-SE,
               names_to = "SF",
               values_to = "weight")

nodes <- data.frame(name = c(tbl$SE, tbl$SF),
                    type = rep(c("SE", "SF"),
                               times = c(length(tbl$SE), length(tbl$SF))))
edges <- data.frame(from = tbl$SF[tbl$weight != 0],
                    to = tbl$SE[tbl$weight != 0])

graph <- tbl_graph(nodes = nodes, edges = edges)



ggraph(graph, 'circlepack') + 
  geom_edge_link() + 
  geom_node_point(aes(colour = type)) +
  coord_fixed()



#~ Plot network ----
#remove unconnected
adj2 <- adj[,colSums(abs(adj)) > 0][1:50,1:50]
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
  wb_tx2g(tx2gt) |> i2s(gids) |> unique()

colSums(abs(adj)) |>
  sort(decreasing = TRUE) |>
  head(n = 22) |>
  names() |>
  wb_tx2g(tx2gt) |>
  i2s(gids) |>
  unique() |>
  clipr::write_clip()


adj_tbl |>
  summarize(`weighted degree` = sum(abs(degree)),
            outdegree = sum(degree != 0),
            .by = c(sf_id, sf_name)) |>
  arrange(desc(`weighted degree`)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = outdegree, y = `weighted degree`))

# outdegree (nb of connections) and weighted degree (sum of connection weights) basically equivalent
adj_tbl |>
  summarize(outdegree = sum(degree != 0),
            .by = c(sf_id, sf_name)) |>
  arrange(desc(outdegree)) |>
  head(20) |>
  clipr::write_clip()



colSums(abs(adj)) |>
  enframe(name = "tx_id", value = "sum_adj") |>
  mutate(gene_id = wb_tx2g(tx_id, tx2gt),
         gene_name = i2s(gene_id, gids)) |>
  summarize(`weighted degree` = sum(sum_adj),
            .by = c(gene_id, gene_name)) |>
  arrange(desc(`weighted degree`)) |>
  head(20) |>
  clipr::write_clip()




# examples ----

#~ cgds-1 ----
sfs_cdgs1 <- adj_tbl_gene |>
  filter(target_name == "cdgs-1") |>
  filter(degree != 0) |>
  pull(sf_id)


adj_tbl_gene |>
  filter(target_name == "cdgs-1") |>
  filter(sf_name %in% c("unc-75"))

se_coords |> filter(gene_id == s2i("cdgs-1", gids))

adj_cdgs1 <- adj["SE_173", wb_g2tx(sfs_cdgs1, tx2gt) |> unlist(), drop = FALSE]
head(adj_cdgs1)




nb_sf <- ncol(adj_cdgs1)
nb_se <- nrow(adj_cdgs1)

adjm <- rbind(
  cbind(matrix(0, nrow = nb_sf, ncol = nb_sf), t(adj_cdgs1)),
  cbind(adj_cdgs1, matrix(0, nrow = nb_se, ncol = nb_se))
)
adjm <- adjm != 0
dimnames(adjm) <- list(c(colnames(adj_cdgs1), rownames(adj_cdgs1)),
                       c(colnames(adj_cdgs1), rownames(adj_cdgs1)))

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
     layout = layout.graphopt(gr, spring.length = 100, spring.constant = .01))





#~ unc-16 ----

adj_tbl_gene |>
  filter(target_name == "unc-16") |>
  filter(degree != 0)

ground_truth_sf <- c("unc-75","exc-7","prp-40")


adj_tbl_gene |>
  filter(target_name == "unc-16",
         degree !=0 |
           sf_name %in% ground_truth_sf) 

# exon coords
# se_coords |>
#   filter(gene_id == s2i("unc-16", gids),
#          startsWith(event_id, "SE")) |>
#   filter(event_id %in% ev_unc16)

sfs_unc16 <- adj_tbl_gene |>
  filter(target_name == "unc-16",
         degree !=0 |
           sf_name %in% ground_truth_sf) |>
  pull(sf_id)

tx_sf_unc16 <- wb_g2tx(sfs_unc16, tx2gt) |> unlist()

adj_tbl_gene |>
  filter(target_name == "cdgs-1") |>
  filter(sf_name %in% c("unc-75"))

ev_unc16 <- se_coords |>
  filter(gene_id == s2i("unc-16", gids),
         startsWith(event_id, "SE")) |>
  pull(event_id) |>
  intersect(rownames(adj))

adj_unc16 <- adj[ev_unc16, wb_g2tx(sfs_unc16, tx2gt) |> unlist(), drop = FALSE]
head(adj_unc16)




nb_sf <- ncol(adj_unc16)
nb_se <- nrow(adj_unc16)

adjm <- rbind(
  cbind(matrix(0, nrow = nb_sf, ncol = nb_sf), t(adj_unc16)),
  cbind(adj_unc16, matrix(0, nrow = nb_se, ncol = nb_se))
)
adjm <- adjm != 0
dimnames(adjm) <- list(c(colnames(adj_unc16), rownames(adj_unc16)),
                       c(colnames(adj_unc16), rownames(adj_unc16)))

gr <- graph_from_adjacency_matrix(adjm,
                                  mode = "directed")



V(gr)$type <- rep(c("SF", "SE"), times = c(nb_sf, nb_se))
V(gr)$shape <- rep(c("circle", "square"), times = c(nb_sf, nb_se))

is_in_ground_truth <- V(gr)[1:nb_sf] |>
  names() |>
  wb_tx2g(tx2gt) |>
  (\(g_id) g_id %in% s2i(ground_truth_sf, gids))()

V(gr)$color <- (c("red4","green3")[is_in_ground_truth+1] |>
  c(rep("blue4", times = nb_se)))



plot(gr,
     vertex.label = NA,
     vertex.size = 4,
     edge.width = 0.1,
     edge.arrow.mode = 0,
     layout = layout.graphopt(gr, spring.length = 100, spring.constant = .01))




#~~ by gene ----
colgenes <- colnames(adj) |> wb_tx2g(tx2gt)

adj_by_gene <- matrix(nrow = nrow(adj),
                      ncol = length(unique(colgenes)),
                      dimnames = list(rownames(adj), unique(colgenes)))
for(cur_gene in unique(colgenes)){
  adj_by_gene[,cur_gene] <- apply(adj[,colgenes == cur_gene,drop=FALSE], 1, biggest)
}
head(adj_by_gene)
head(adj)

colnames(adj_by_gene) <- i2s(colnames(adj_by_gene), gids, warn_missing = TRUE)






ground_truth_sf <- c("unc-75","exc-7","prp-40")


adj_tbl_gene |>
  filter(target_name == "unc-16",
         degree !=0 |
           sf_name %in% ground_truth_sf) 




sfs_unc16 <- adj_tbl_gene |>
  filter(target_name == "unc-16",
         degree !=0 |
           sf_name %in% ground_truth_sf) |>
  pull(sf_id)

sf_unc16 <- i2s(sfs_unc16, gids)



ev_unc16 <- se_coords |>
  filter(gene_id == s2i("unc-16", gids),
         startsWith(event_id, "SE")) |>
  pull(event_id) |>
  intersect(rownames(adj))



adj_unc16 <- adj_by_gene[ev_unc16, sf_unc16]
head(adj_unc16)




nb_sf <- ncol(adj_unc16)
nb_se <- nrow(adj_unc16)

adjm <- rbind(
  cbind(matrix(0, nrow = nb_sf, ncol = nb_sf), t(adj_unc16)),
  cbind(adj_unc16, matrix(0, nrow = nb_se, ncol = nb_se))
)
adjm <- adjm != 0
dimnames(adjm) <- list(c(colnames(adj_unc16), rownames(adj_unc16)),
                       c(colnames(adj_unc16), rownames(adj_unc16)))

gr <- graph_from_adjacency_matrix(adjm,
                                  mode = "undirected")



V(gr)$type <- rep(c("SF", "SE"), times = c(nb_sf, nb_se))
V(gr)$shape <- rep(c("circle", "square"), times = c(nb_sf, nb_se))

is_in_ground_truth <- V(gr)[1:nb_sf] |>
  names() |>
  (\(g_id) g_id %in% ground_truth_sf)()

V(gr)$color <- (c("red4","green3")[is_in_ground_truth+1] |>
                  c(rep("blue4", times = nb_se)))



plot(gr,
     vertex.label = NA,
     vertex.size = 4,
     edge.width = 0.1,
     edge.arrow.mode = 0)



RCy3::createNetworkFromIgraph(gr,"unc16")


# --> source data
# trgt <- "unc-16"
# se_coords |>
#   filter(gene_id == s2i(trgt, gids),
#          startsWith(event_id, "SE")) |>
#   clipr::write_clip()
# 
# adj_tbl |>
#   filter(target_name == trgt,
#          degree !=0 |
#            sf_name %in% ground_truth_sf) |>
#   clipr::write_clip()





#~ lin-10 ----

adj_tbl_gene |>
  filter(target_name == "lin-10") |>
  filter(degree != 0)


#~ C07A12.7 ----

adj_tbl_gene |>
  filter(target_name == "C07A12.7") |>
  filter(degree != 0) |>
  pull(sf_name) |> sort() |> paste0(collapse = ", ") |> writeClipboard()



#~~ by gene ----
trgt <- "C07A12.7"
ground_truth_sf <- c("unc-75")


adj_tbl_gene |>
  filter(target_name == trgt,
         degree !=0 |
           sf_name %in% ground_truth_sf) 




sfs_sub <- adj_tbl_gene |>
  filter(target_name == trgt,
         degree !=0 |
           sf_name %in% ground_truth_sf) |>
  pull(sf_id)

sfs_sub <- i2s(sfs_sub, gids)



ev_sub <- se_coords |>
  filter(gene_id == s2i(trgt, gids),
         startsWith(event_id, "SE")) |>
  pull(event_id) |>
  intersect(rownames(adj))



adj_sub <- adj_by_gene[ev_sub, sfs_sub, drop = FALSE]
head(adj_sub)




nb_sf <- ncol(adj_sub)
nb_se <- nrow(adj_sub)

adjm <- rbind(
  cbind(matrix(0, nrow = nb_sf, ncol = nb_sf), t(adj_sub)),
  cbind(adj_sub, matrix(0, nrow = nb_se, ncol = nb_se))
)
adjm <- adjm != 0
dimnames(adjm) <- list(c(colnames(adj_sub), rownames(adj_sub)),
                       c(colnames(adj_sub), rownames(adj_sub)))

gr <- graph_from_adjacency_matrix(adjm,
                                  mode = "undirected")



V(gr)$type <- rep(c("SF", "SE"), times = c(nb_sf, nb_se))
V(gr)$shape <- rep(c("circle", "square"), times = c(nb_sf, nb_se))

is_in_ground_truth <- V(gr)[1:nb_sf] |>
  names() |>
  (\(g_id) g_id %in% ground_truth_sf)()

V(gr)$color <- (c("red4","green3")[is_in_ground_truth+1] |>
                  c(rep("blue4", times = nb_se)))



plot(gr,
     vertex.label = NA,
     vertex.size = 10,
     edge.width = 0.1,
     edge.arrow.mode = 0)


RCy3::cytoscapePing()
RCy3::createNetworkFromIgraph(gr,"C07A12.7")



# --> source data
# se_coords |>
#   filter(gene_id == s2i(trgt, gids),
#          startsWith(event_id, "SE")) |>
#   clipr::write_clip()
# 
# adj_tbl |>
#   filter(target_name == trgt,
#          degree !=0 |
#            sf_name %in% ground_truth_sf) |>
#   clipr::write_clip()


#~ ret-1 -----

adj_tbl_gene |>
  filter(target_name == "ret-1") |>
  filter(degree != 0)

#~~ by gene ----
trgt <- "ret-1"
ground_truth_sf <- c("unc-75", "hrpr-1","sfa-1","uaf-2","snr-1","prp-38")


adj_tbl_gene |>
  filter(target_name == trgt,
         degree !=0 |
           sf_name %in% ground_truth_sf) 




sfs_sub <- adj_tbl_gene |>
  filter(target_name == trgt,
         degree !=0 |
           sf_name %in% ground_truth_sf) |>
  pull(sf_id)

sfs_sub <- i2s(sfs_sub, gids)



ev_sub <- se_coords |>
  filter(gene_id == s2i(trgt, gids),
         startsWith(event_id, "SE")) |>
  pull(event_id) |>
  intersect(rownames(adj))



adj_sub <- adj_by_gene[ev_sub, sfs_sub, drop = FALSE]
head(adj_sub)




nb_sf <- ncol(adj_sub)
nb_se <- nrow(adj_sub)

adjm <- rbind(
  cbind(matrix(0, nrow = nb_sf, ncol = nb_sf), t(adj_sub)),
  cbind(adj_sub, matrix(0, nrow = nb_se, ncol = nb_se))
)
adjm <- adjm != 0
dimnames(adjm) <- list(c(colnames(adj_sub), rownames(adj_sub)),
                       c(colnames(adj_sub), rownames(adj_sub)))

gr <- graph_from_adjacency_matrix(adjm,
                                  mode = "undirected")



V(gr)$type <- rep(c("SF", "SE"), times = c(nb_sf, nb_se))
V(gr)$shape <- rep(c("circle", "square"), times = c(nb_sf, nb_se))

is_in_ground_truth <- V(gr)[1:nb_sf] |>
  names() |>
  (\(g_id) g_id %in% ground_truth_sf)()

V(gr)$color <- (c("red4","green3")[is_in_ground_truth+1] |>
                  c(rep("blue4", times = nb_se)))



plot(gr,
     vertex.label = NA,
     vertex.size = 4,
     edge.width = 0.1,
     edge.arrow.mode = 0)


RCy3::cytoscapePing()
RCy3::createNetworkFromIgraph(gr,"re-1")




#~ daf-2 ----


adj_tbl_gene |>
  filter(target_name == "daf-2") |>
  filter(degree != 0) |>
  pull(sf_name) |>
  setdiff(ground_truth_sf) |>
  sort() |>
  paste0(collapse = ", ") |>
  writeClipboard()

adj_tbl_gene |>
  filter(target_name == "daf-2") |>
  filter(sf_name %in% c("unc-75", "ptb-1", "asd-1", "hrpa-1", "hrpf-1", "exc-7", "rsp-8"))





#~~ by gene ----
trgt <- "daf-2"
ground_truth_sf <- c("unc-75", "ptb-1", "asd-1", "hrpa-1",
                     "hrpf-1", "exc-7", "rsp-8", "rsp-2")


adj_tbl_gene |>
  filter(target_name == trgt,
         degree !=0 |
           sf_name %in% ground_truth_sf) 




sfs_sub <- adj_tbl_gene |>
  filter(target_name == trgt,
         degree !=0 |
           sf_name %in% ground_truth_sf) |>
  pull(sf_id)

sfs_sub <- i2s(sfs_sub, gids)



ev_sub <- se_coords |>
  filter(gene_id == s2i(trgt, gids),
         startsWith(event_id, "SE")) |>
  pull(event_id) |>
  intersect(rownames(adj))



adj_sub <- adj_by_gene[ev_sub, sfs_sub, drop = FALSE]
head(adj_sub)




nb_sf <- ncol(adj_sub)
nb_se <- nrow(adj_sub)

adjm <- rbind(
  cbind(matrix(0, nrow = nb_sf, ncol = nb_sf), t(adj_sub)),
  cbind(adj_sub, matrix(0, nrow = nb_se, ncol = nb_se))
)
adjm <- adjm != 0
dimnames(adjm) <- list(c(colnames(adj_sub), rownames(adj_sub)),
                       c(colnames(adj_sub), rownames(adj_sub)))

gr <- graph_from_adjacency_matrix(adjm,
                                  mode = "undirected")



V(gr)$type <- rep(c("SF", "SE"), times = c(nb_sf, nb_se))
V(gr)$shape <- rep(c("circle", "square"), times = c(nb_sf, nb_se))

is_in_ground_truth <- V(gr)[1:nb_sf] |>
  names() |>
  (\(g_id) g_id %in% ground_truth_sf)()

V(gr)$color <- (c("red4","green3")[is_in_ground_truth+1] |>
                  c(rep("blue4", times = nb_se)))



plot(gr,
     vertex.label = NA,
     vertex.size = 4,
     edge.width = 0.1,
     edge.arrow.mode = 0)


RCy3::cytoscapePing()
RCy3::createNetworkFromIgraph(gr,"daf-2")



# --> source data
se_coords |>
  filter(gene_id == s2i(trgt, gids),
         startsWith(event_id, "SE")) |>
  clipr::write_clip()

adj_tbl |>
  filter(target_name == trgt,
         degree !=0 |
           sf_name %in% ground_truth_sf) |>
  clipr::write_clip()


