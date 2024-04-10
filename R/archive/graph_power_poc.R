# Predictive power of graph approach
#
# Using out-of-sample: use some of the samples to make the graph, use the rest of the samples to predict PSI
# Use DPM for the graph



library(tidyverse)
library(DPM)
library(wbData)

gids <- wb_load_gene_ids(281)

quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs") |>
  filter(transcript_id != "R07E5.14.2")

# events_coords <- read_tsv("data/export_for_arman/221111_events_coordinates.tsv") |>
#   select(event_id, gene_id)
# 
# list_sfs <- unique(sf_expression$gene_id)



# Prepare data ----

#~ PSI -----
mat_psi <- quantifs_filtered |>
  mutate(dPSI = PSI - mean(PSI),
         .by = event_id) |>
  select(event_id, sample_id, dPSI) |>
  pivot_wider(id_cols = sample_id, names_from = event_id, values_from = dPSI) |>
  column_to_rownames("sample_id") |>
  as.matrix()

# filter PSI
dim(mat_psi)
mat_psi[1:3,1:3]


# remove samples full of NA
rowMeans(is.na(mat_psi)) |> hist(breaks = 50)
mat_psi <- mat_psi[rowMeans(is.na(mat_psi)) < .4, ]

dim(mat_psi)
mat_psi[1:3,1:3]

# Train/test split: note we do that BEFORE imputation
set.seed(123)
train_samples <- sample(rownames(mat_psi), size = .7*nrow(mat_psi))
test_samples <- setdiff(rownames(mat_psi), train_samples)

mat_psi_train <- mat_psi[train_samples,]


# Impute the rest of the NAs
# note: Park, Wang and Lim (2020) https://arxiv.org/abs/2006.04632 try a few imputation methods, {impute} seems OK
mat_psi_train <- impute::impute.knn(mat_psi_train)$data
dim(mat_psi_train)
mat_psi_train[1:3,1:3]



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


dim(mat_psi_train)
dim(mat_sf_train)
# match rows
all.equal(rownames(mat_psi_train), rownames(mat_sf_train))

mat_psi_train[1:3,1:3]
mat_sf_train[1:3,1:3]


mat_train <- cbind(mat_sf_train, mat_psi_train)


dim(mat_train)
mat_train[1:3,1:3]

# finish
mat_psi_test <- mat_psi[test_samples, ]
mat_sf_test <- mat_sf[test_samples,]
nb_psi <- ncol(mat_psi_train)
nb_sf <- ncol(mat_sf_train)




# Make network ----

dpm <- DPM::reg.dpm(mat_train)

#  adjacency matrix
adj <- dpm[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]

rownames(adj) <- colnames(mat_psi_train)
colnames(adj) <- colnames(mat_sf_train)

dim(adj)
adj[1:3,1:3]

image(adj)
pheatmap::pheatmap(adj)


# Evaluate on test set ----

dim(mat_sf_test)
mat_sf_test[1:3,1:3]

predicted_psi <- mat_sf_test %*% t(adj)

dim(predicted_psi)
predicted_psi[1:3,1:3]
mat_psi_test[1:3,1:3]

plot(mat_psi_test, predicted_psi)

plot(mat_psi_test, abs(predicted_psi - mat_psi_test))


# Estimate deviance
rmse <- sqrt(sum(predicted_psi^2 - mat_psi_test^2, na.rm = TRUE)/length(mat_psi_test))
mae <- sum(abs(predicted_psi - mat_psi_test), na.rm = TRUE)/length(mat_psi_test)





# Non-regularized DPM ----

dpm <- DPM::reg.dpm(mat_train)

#  adjacency matrix
adj <- dpm[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]

predicted_psi <- mat_sf_test %*% t(adj)

mae <- sum(abs(predicted_psi - mat_psi_test), na.rm = TRUE)/length(mat_psi_test)
mae

plot(mat_psi_test, predicted_psi)
plot(mat_psi_test, abs(predicted_psi - mat_psi_test))



# ARACNE ----

mim <- minet::build.mim(mat_train,
                 estimator = "spearman",
                 disc = "none")
arac <- minet::aracne(mim, eps = 0)

adj <- arac[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]

predicted_psi <- mat_sf_test %*% t(adj)

mae <- sum(abs(predicted_psi - mat_psi_test), na.rm = TRUE)/length(mat_psi_test)
mae

pheatmap::pheatmap(adj, show_rownames = FALSE, show_colnames = FALSE)
plot(mat_psi_test, predicted_psi)
plot(mat_psi_test, abs(predicted_psi - mat_psi_test))

# opar <- par()
par(mfrow = c(3,1), mar = c(3,0.5,0,.5))
j <- sample(rownames(mat_psi_test), 1)
plot(mat_psi_test[j,], predicted_psi[j,])
i <- sample(colnames(mat_psi_test), 1)
plot(mat_psi_test[,i], predicted_psi[,i])


# CV on eps
get_arac <- function(eps){
  arac <- minet::aracne(mim, eps)
  
  adj <- arac[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]
  
  predicted_psi <- mat_sf_test %*% t(adj)
  
  mae <- sum(abs(predicted_psi - mat_psi_test), na.rm = TRUE)/(length(mat_psi_test) - sum(is.na(mat_psi_test)))
  mae
}


tibble(eps = (0:10)/50,
       mae = map_dbl(eps, get_arac)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = eps, y = mae))




# more ARACNE

mim <- minet::build.mim(mat_train,
                        estimator = "pearson",
                        disc = "none")
arac <- minet::aracne(mim, eps = 0)

adj <- arac[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]

predicted_psi <- mat_sf_test %*% t(adj)

mae <- sum(abs(predicted_psi - mat_psi_test), na.rm = TRUE)/length(mat_psi_test)
mae

pheatmap::pheatmap(adj, show_rownames = FALSE, show_colnames = FALSE)
plot(mat_psi_test, predicted_psi)
plot(mat_psi_test, abs(predicted_psi - mat_psi_test))


# more more ARACNE

mim <- minet::build.mim(mat_train,
                        estimator = "mi.empirical",
                        disc = "equalfreq")
arac <- minet::aracne(mim, eps = 0)

adj <- arac[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]

predicted_psi <- mat_sf_test %*% t(adj)

mae <- sum(abs(predicted_psi - mat_psi_test), na.rm = TRUE)/length(mat_psi_test)
mae

pheatmap::pheatmap(adj, show_rownames = FALSE, show_colnames = FALSE)
plot(mat_psi_test, predicted_psi)
plot(mat_psi_test, abs(predicted_psi - mat_psi_test))


# more mrnet

mim <- minet::build.mim(mat_train,
                        estimator = "spearman",
                        disc = "none")
arac <- minet::mrnet(mim)

adj <- arac[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]

predicted_psi <- mat_sf_test %*% t(adj)

mae <- sum(abs(predicted_psi - mat_psi_test), na.rm = TRUE)/length(mat_psi_test)
mae

pheatmap::pheatmap(adj, show_rownames = FALSE, show_colnames = FALSE)
plot(mat_psi_test, predicted_psi)
plot(mat_psi_test, abs(predicted_psi - mat_psi_test))










