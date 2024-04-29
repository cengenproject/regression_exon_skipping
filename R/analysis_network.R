# Network analysis
#
# From cluster

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
  filter(sf_name == "asd-1") #|> filter(degree != 0)
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


# by target
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

test_res



#~ SF GO ----

adj_tbl_gene$sf_id |>
  unique() |>
  writeLines("presentations/240308_figures/background.txt")

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
  wb_tx2g(tx2gt) |> i2s(gids)


