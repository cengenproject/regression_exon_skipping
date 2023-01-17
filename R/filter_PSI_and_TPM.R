# Filtering and cleaning for lasso on PSI vs TPM

# Packages ----

library(tidyverse)
# library(wbData)


# tx2g <- wb_load_tx2gene(281)
# gids <- wb_load_gene_ids(281)

sample2neuron <- function(sample_ids){
  res <- str_extract(sample_ids, "^([A-Z0-9ef]{2,4})r\\d{1,3}$", group = 1L)
  if(any(is.na(res))) warning(sum(is.na(res)), " NAs returned: ", head(sample_ids[is.na(res)]))
  res
}


# Read data ----

quantifs <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Z1-9]{2,4})r[0-9]{2,3}$")[,2]) 

all_tx_expression <- read_tsv("data/export_for_arman/tx_expression.tsv.gz") 




# Filter neurons ----

# Filter neurons with too few samples. Also remove Ref
keep_neurons <- all_tx_expression |>
  select(sample_id, neuron_id) |>
  distinct() |>
  dplyr::count(neuron_id) |>
  filter(n > 2) |>
  pull(neuron_id) |>
  setdiff("Ref")


all_tx_expression <- all_tx_expression |>
  filter(neuron_id %in% keep_neurons)

quantifs <- quantifs |>
  filter(neuron_id %in% keep_neurons)


# Re-Make tx matrix

mat_tx_expr <- all_tx_expression |>
  pivot_wider(transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr <- mat_tx_expr[rowSums(mat_tx_expr) != 0,]
dim(mat_tx_expr)





#~ Filter samples based on tx data ----

# PCA
pc <- prcomp(t(log1p(mat_tx_expr)), center = TRUE, scale. = TRUE)


pcdf <- cbind(as.data.frame(pc$x), data.frame(sample_id = colnames(mat_tx_expr))) |>
  mutate(neuron_id = sample2neuron(sample_id))

plot(pcdf[,1:2])
ggplot(pcdf, aes(x = PC1 , y = PC2, color = neuron_id)) +
  theme_classic() +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id)) +
  geom_abline(aes(intercept = -150, slope = -.3)) +
  hues::scale_color_iwanthue()

keep_samples <- unique(all_tx_expression$sample_id) |>
  setdiff(c("RICr133","PVMr122","AVKr113"))


# Filter bad samples
all_tx_expression <- all_tx_expression |>
  filter(sample_id %in% keep_samples)

quantifs <- quantifs |>
  filter(sample_id %in% keep_samples)


# Re-Make tx matrix
mat_tx_expr <- all_tx_expression |>
  pivot_wider(transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr <- mat_tx_expr[rowSums(mat_tx_expr) != 0,]
dim(mat_tx_expr)






#~ Filter samples based on exon data ----

# make matrix
mat_quantifs <- quantifs |>
  pivot_wider(event_id,
              names_from = "sample_id",
              values_from = "PSI") |>
  column_to_rownames("event_id") |>
  as.matrix()
dim(mat_quantifs)

# need to tackle NAs before PCA
hist(rowSums(is.na(mat_quantifs)))
hist(colSums(is.na(mat_quantifs)))

mat_quantifs <- mat_quantifs[rowSums(is.na(mat_quantifs)) < 100,
                             colSums(is.na(mat_quantifs)) < 500]
dim(mat_quantifs)

mat_quantifs <- missMDA::imputePCA(mat_quantifs)
dim(mat_quantifs$completeObs)

# PCA
pc <- prcomp(t(mat_quantifs$completeObs), center = TRUE, scale. = TRUE)


pcdf <- cbind(as.data.frame(pc$x), data.frame(sample_id = colnames(mat_quantifs$completeObs))) |>
  mutate(neuron_id = sample2neuron(sample_id))

plot(pcdf[,1:2])
ggplot(pcdf, aes(x = PC1 , y = PC2, color = neuron_id)) +
  theme_classic() +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id)) +
  hues::scale_color_iwanthue()


#> not filtering anything else












# Filter tx ----

# PCA on tx
pc <- prcomp(log1p(mat_tx_expr), center = TRUE, scale. = TRUE)


pcdf <- cbind(as.data.frame(pc$x), data.frame(tx_id = rownames(mat_tx_expr)))
plot(pcdf[,1:2])
abline(v = 80)
pcdf[pcdf$PC1 > 60,"tx_id"] |> writeClipboard()
#> ribosomal, mitochondrial, ncRNAs, etc


# table(wbData::wb_tx2g(pcdf[pcdf$PC1 > 60,"tx_id"], tx2g) %in% wormDatasets::worm_putative_splice_factors$gene_id)

keep_tx <- rownames(mat_tx_expr) |>
  setdiff(pcdf[pcdf$PC1 > 60,"tx_id"])



# Filter bad tx
all_tx_expression <- all_tx_expression |>
  filter(transcript_id %in% keep_tx)


# Re-Make tx matrix
mat_tx_expr <- all_tx_expression |>
  pivot_wider(transcript_id,
              names_from = "sample_id",
              values_from = "TPM") |>
  column_to_rownames("transcript_id") |>
  as.matrix()

# remove empty rows
mat_tx_expr <- mat_tx_expr[rowSums(mat_tx_expr) != 0,]
dim(mat_tx_expr)







#~ Filter exons ----

# make matrix
mat_quantifs <- quantifs |>
  pivot_wider(event_id,
              names_from = "sample_id",
              values_from = "PSI") |>
  column_to_rownames("event_id") |>
  as.matrix()
dim(mat_quantifs)

# need to tackle NAs before PCA
hist(rowSums(is.na(mat_quantifs)))
hist(colSums(is.na(mat_quantifs)))

mat_quantifs <- mat_quantifs[rowSums(is.na(mat_quantifs)) < 100,
                             colSums(is.na(mat_quantifs)) < 500]
dim(mat_quantifs)

mat_quantifs <- missMDA::imputePCA(mat_quantifs)
dim(mat_quantifs$completeObs)

# PCA
pc <- prcomp(mat_quantifs$completeObs, center = TRUE, scale. = TRUE)


pcdf <- cbind(as.data.frame(pc$x), data.frame(sample_id = rownames(mat_quantifs$completeObs)))

plot(pcdf[,1:2])


#> not filtering anything else






# Save filtered values ----



# qs::qsave(quantifs, "data/intermediates/230117_quantifs_filtered.qs")
# 
# qs::qsave(all_tx_expression, "data/intermediates/230117_tx_filtered.qs")





