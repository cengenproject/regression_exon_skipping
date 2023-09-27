# Predictive power of graph approach
#
# Using out-of-sample: use some of the samples to make the graph, use the rest of the samples to predict PSI
# Use DPM for the graph. After the Proof of Concept script, use a more systematic cross-validation


# Inits ----

library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)

source("R/loss_functions.R")



quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs") |>
  filter(transcript_id != "R07E5.14.2")


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


# remove samples full of NA
mat_psi <- mat_psi[rowMeans(is.na(mat_psi)) < .4, ]


# Train/test split: note we do that BEFORE imputation
set.seed(123)
train_samples <- sample(rownames(mat_psi), size = .7*nrow(mat_psi))
test_samples <- setdiff(rownames(mat_psi), train_samples)

mat_psi_train <- mat_psi[train_samples,]


# Impute the rest of the NAs
# note: Park, Wang and Lim (2020) https://arxiv.org/abs/2006.04632 try a few imputation methods, {impute} seems OK
mat_psi_train <- impute::impute.knn(mat_psi_train)$data


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



# match rows
stopifnot(all.equal(rownames(mat_psi_train), rownames(mat_sf_train)))
mat_train <- cbind(mat_sf_train, mat_psi_train)


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









# Simulated Precision matrix ----
# First question: which do act as inverse of cov? when cov invertible
set.seed(1)
gen <- huge::huge.generator(graph = "scale-free",vis = TRUE); par(mfrow = c(1,1))
Y <- gen$data

S <- cov(Y)

image(Y)
image(S)

#~ Direct inversion ----
OM <- solve(S)
image(OM)
S %*% OM
image(S %*% OM)

#~ DPM ----
OM <- DPM::reg.dpm(Y)
S %*% OM
image(S %*% OM)


#~ QUIC ----

OM <- QUIC::QUIC(S, rho = matrix(1, nrow = nrow(S), ncol = ncol(S)))

S %*% OM$X
image(S %*% OM$X)


#~ glasso ----
OM <- glasso::glasso(S, rho = 1)

S %*% OM$wi
image(S %*% OM$wi)


#~ sparseMatEst ----
OM <- sparseMatEst::sparsePrec(Y)
image(S %*% OM$hard[,,3])

#~ clime ----
library(clime)
OM <- clime::clime(Y)
OM2 <- clime::cv.clime(OM)
OM <- clime(Y, lambda = OM2$lambdaopt)

image(S %*% OM$Omegalist[[1]])


#~ scio ----
OM <- scio::scio(S, 0.5)
image(S %*% OM$w)


#~ huge ----
OM <- huge::huge.npn(Y) |>
  huge::huge() |>
  huge::huge.select()

plot(OM)
huge:::plot.huge(OM)



# Compute loss functions ----



#~ glasso ----


cv_glasso_res <- map(c(2, 1, .5, .125, .1,.08, .075, .06, .05, .04, .01, .005, .001) |> set_names(),
                     \(.rho){
                       message(.rho)
                       map(unique(folds) |> set_names(),
                           \(fold){
                             S_train <- cov(mat_train[folds != fold,])
                             S_valid <- cov(mat_train[folds == fold,])
                             
                             OM_train <- glasso::glasso(S_train, rho = .rho)
                             
                             list(S_train = S_train,
                                  S_valid = S_valid,
                                  res = OM_train)
                           },
                           .progress = TRUE)
                     })

# qs::qsave(cv_glasso_res, "data/intermediates/230818_cv/230824_cv_glasso_res.qs")
cv_glasso_res <- qs::qread("data/intermediates/230818_cv/230824_cv_glasso_res.qs")

cv_glasso <- cv_glasso_res |>
  setNames(paste0(names(cv_glasso_res), "|")) |>
  unlist(recursive = FALSE) |>
  enframe() |>
  separate_wider_delim(name,
                       delim = "|.",
                       names = c("rho", "fold")) |>
  unnest_wider(value) |>
  mutate(rho = as.numeric(rho),
         fold = as.integer(fold))


cv_glasso <- cv_glasso |>
  mutate(loss_frob = map2_dbl(S_valid, res,
                              ~ loss_frob(.x, .y[["w"]])),
         loss_quad = map2_dbl(S_valid, res,
                              ~ loss_quad(.x, .y[["wi"]])),
         loss_cv_I = map2_dbl(fold, res,
                              ~loss_tong_cv_I(mat_train[folds != .x,], .y[["wi"]])),
         loss_cv_II = map2_dbl(fold, res,
                               ~loss_tong_cv_II(mat_train[folds != .x,], .y[["wi"]]))) |>
  mutate(adj_sparsity = map_dbl(res,
                                ~ mat_sparsity(.x[["wi"]])),
         adj_power_law_R2 = map_dbl(res,
                                    ~ mat_power_law(.x[["wi"]])))


# qs::qsave(cv_glasso, "data/intermediates/230818_cv/230828_cv_glasso.qs")
cv_glasso <- qs::qread("data/intermediates/230818_cv/230828_cv_glasso.qs")



cv_glasso |>
  select(where(is.atomic)) |>
  pivot_longer(cols = starts_with("loss"),
               names_prefix = "loss_",
               names_to = "loss_function",
               values_to = "loss") |>
  group_by(rho, loss_function) |>
  summarize(mean_loss = mean(loss),
            sd_loss = sd(loss),
            .groups = 'drop') |>
  ggplot(aes(x = rho, y = mean_loss)) +
  theme_classic() +
  facet_grid(rows=vars(loss_function), scales = "free_y") +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_loss - sd_loss, ymax = mean_loss + sd_loss),
                width = .05) +
  scale_x_log10() 

cv_glasso |>
  select(where(is.atomic)) |>
  pivot_longer(cols = starts_with("adj"),
               names_prefix = "adj_",
               names_to = "loss_function",
               values_to = "metric") |>
  group_by(rho, loss_function) |>
  summarize(mean_metric = mean(metric),
            sd_metric = sd(metric),
            .groups = 'drop') |>
  ggplot(aes(x = rho, y = mean_metric)) +
  theme_classic() +
  facet_grid(rows=vars(loss_function), scales = "free_y") +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_metric - sd_metric, ymax = mean_metric + sd_metric),
                width = .05) +
  geom_hline(aes(yintercept = 1), color = 'grey') +
  scale_x_log10() 






#~ QUIC ----

run_quic <- function(S, r_gg, r_gs, r_ss){
  regul_mat <- rbind(cbind(matrix(r_gg, nrow=nb_sf, ncol = nb_sf),
                           matrix(r_gs, nrow=nb_sf, ncol = nb_psi)),
                     cbind(matrix(r_gs, nrow=nb_psi, ncol = nb_sf),
                           matrix(r_ss, nrow=nb_psi, ncol = nb_psi)))
  
  QUIC::QUIC(S, rho = regul_mat, msg = 0)
}


cv_quic_res <- tibble(r_gg = c(.01, .025, .05, .1, .5, 1, 2),
                      r_gs = r_gg,
                      r_ss = r_gg) |>
  pmap(\(r_gg, r_gs, r_ss){
    
    message(r_gg,"; ", r_gs,"; ", r_ss)
    map(unique(folds) |> set_names(),
        \(fold){
          
          S_train <- cov(mat_train[folds != fold,])
          S_valid <- cov(mat_train[folds == fold,])
          
          OM_train <- run_quic(S_train, r_gg, r_gs, r_ss)
          
          list(S_train = S_train,
               S_valid = S_valid,
               res = OM_train)
        },
        .progress = TRUE)
  })

# qs::qsave(cv_quic_res, "data/intermediates/230818_cv/230824_cv_quic_res.qs")
cv_quic_res <- qs::qread("data/intermediates/230818_cv/230824_cv_quic_res.qs")

cv_quic <- cv_quic_res |>
  unlist(recursive = FALSE) |>
  enframe() |>
  bind_cols(expand_grid(rho = c(.01, .025, .05, .1, .5, 1, 2),
                        fold = unique(folds))) |>
  unnest_wider(value) |>
  mutate(loss_frob = map2_dbl(S_valid, res,
                              ~ loss_frob(.x, .y[["W"]])),
         loss_quad = map2_dbl(S_valid, res,
                              ~ loss_quad(.x, .y[["X"]])),
         loss_cv_I = map2_dbl(fold, res,
                              ~loss_tong_cv_I(mat_train[folds != .x,], .y[["X"]])),
         loss_cv_II = map2_dbl(fold, res,
                               ~loss_tong_cv_II(mat_train[folds != .x,], .y[["X"]]))) |>
  mutate(adj_sparsity = map_dbl(res,
                            ~ mat_sparsity(.x[["X"]])),
         adj_power_law_R2 = map_dbl(res,
                                ~ mat_power_law(.x[["X"]])))

# qs::qsave(cv_quic, "data/intermediates/230818_cv/230828_cv_quic.qs")
cv_quic <- qs::qread("data/intermediates/230818_cv/230828_cv_quic.qs")


cv_quic |>
  select(where(is.atomic)) |>
  pivot_longer(cols = starts_with("loss"),
               names_prefix = "loss_",
               names_to = "loss_function",
               values_to = "loss") |>
  group_by(rho, loss_function) |>
  summarize(mean_loss = mean(loss),
            sd_loss = sd(loss),
            .groups = 'drop') |>
  ggplot(aes(x = rho, y = mean_loss)) +
  theme_classic() +
  facet_grid(rows=vars(loss_function), scales = "free_y") +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_loss - sd_loss, ymax = mean_loss + sd_loss),
                width = .05) +
  scale_x_log10() 

cv_quic |>
  select(where(is.atomic)) |>
  pivot_longer(cols = starts_with("adj"),
               names_prefix = "adj_",
               names_to = "loss_function",
               values_to = "metric") |>
  group_by(rho, loss_function) |>
  summarize(mean_metric = mean(metric),
            sd_metric = sd(metric),
            .groups = 'drop') |>
  ggplot(aes(x = rho, y = mean_metric)) +
  theme_classic() +
  facet_grid(rows=vars(loss_function), scales = "free_y") +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_metric - sd_metric, ymax = mean_metric + sd_metric),
                width = .05) +
  geom_hline(aes(yintercept = 1), color = 'grey') +
  scale_x_log10() 




#~ Note: compare QUIC and glasso ----

cv_glasso |>
  filter(rho == 0.05) |>
  select(where(is.atomic)) |>
  summarize(across(-c(rho, fold),
                   mean))
  
cv_quic |>
  filter(rho == 0.05) |>
  select(where(is.atomic)) |>
  summarize(across(-c(rho, fold),
                   mean))


corpcor::is.positive.definite(cv_quic$res[[which(cv_quic$rho == 0.05 & cv_quic$fold == 1)]][["X"]])
corpcor::is.positive.definite(cv_glasso$res[[which(cv_glasso$rho == 0.05 & cv_glasso$fold == 1)]][["wi"]])


#~ scio ----

cv_scio_res <- map(c(100, 50, 10, 5, 1, .5, .1, .05, .01) |> set_names(),
               \(.lambda){
                 message(.lambda)
                 map(unique(folds) |> set_names(),
                     \(fold){
                       S_train <- cov(mat_train[folds != fold,])
                       S_valid <- cov(mat_train[folds == fold,])
                       
                       OM_train <- scio::scio(S_train, lambda = .lambda)
                       
                       list(S_train = S_train,
                            S_valid = S_valid,
                            res = OM_train)
                     },
                     .progress = TRUE)
               })


# qs::qsave(cv_scio_res, "data/intermediates/230818_cv/23082_cv_scio_res.qs")
cv_scio_res <- qs::qread("data/intermediates/230818_cv/23082_cv_scio_res.qs")


xx <- map(cv_scio_res,
    \(.pen) map_lgl(.pen,
          \(.fold) isSymmetric(.fold[["res"]][["w"]])))


cv_scio <- cv_scio_res |>
  setNames(paste0(names(cv_scio_res), "|")) |>
  unlist(recursive = FALSE) |>
  enframe() |>
  separate_wider_delim(name,
                       delim = "|.",
                       names = c("rho", "fold")) |>
  unnest_wider(value) |>
  mutate(rho = as.numeric(rho),
         fold = as.integer(fold)) |>
  mutate(res = map(res, ~corpcor::make.positive.definite(.x[["w"]])))

cv_scio2 <- cv_scio |>
  mutate(loss_frob = map2_dbl(S_valid, res,
                              ~ loss_frob(.x, corpcor::pseudoinverse(.y))),
         loss_quad = map2_dbl(S_valid, res,
                              ~ loss_quad(.x, .y)),
         loss_cv_I = map2_dbl(fold, res,
                              ~loss_tong_cv_I(mat_train[folds != .x,], .y)),
         loss_cv_II = map2_dbl(fold, res,
                               ~loss_tong_cv_II(mat_train[folds != .x,], .y))) |>
  mutate(adj_sparsity = map_dbl(res,
                                ~ mat_sparsity(.x)),
         adj_power_law_R2 = map_dbl(res,
                                    ~ mat_power_law(.x)))


# qs::qsave(cv_scio2, "data/intermediates/230818_cv/230829_cv_scio.qs")
cv_scio2 <- qs::qread("data/intermediates/230818_cv/230829_cv_scio.qs")




cv_scio2 |>
  select(where(is.atomic)) |>
  pivot_longer(cols = starts_with("loss"),
               names_prefix = "loss_",
               names_to = "loss_function",
               values_to = "loss") |>
  group_by(rho, loss_function) |>
  summarize(mean_loss = mean(loss),
            sd_loss = sd(loss),
            .groups = 'drop') |>
  ggplot(aes(x = rho, y = mean_loss)) +
  theme_classic() +
  facet_grid(rows=vars(loss_function), scales = "free_y") +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_loss - sd_loss, ymax = mean_loss + sd_loss),
                width = .05) +
  scale_x_log10()

cv_scio2 |>
  select(where(is.atomic)) |>
  pivot_longer(cols = starts_with("adj"),
               names_prefix = "adj_",
               names_to = "loss_function",
               values_to = "metric") |>
  group_by(rho, loss_function) |>
  summarize(mean_metric = mean(metric, na.rm = TRUE),
            sd_metric = sd(metric, na.rm = TRUE),
            .groups = 'drop') |>
  ggplot(aes(x = rho, y = mean_metric)) +
  theme_classic() +
  facet_grid(rows=vars(loss_function), scales = "free_y") +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_metric - sd_metric, ymax = mean_metric + sd_metric),
                width = .05) +
  geom_hline(aes(yintercept = 1), color = 'grey') +
  scale_x_log10() 


#~ SCIO provided CV by likelihood ----

cv_sciocv_res <- map(unique(folds) |> set_names(),
                   \(fold){
                     S_train <- cov(mat_train[folds != fold,])
                     S_valid <- cov(mat_train[folds == fold,])
                     
                     OM_train <- scio::scio.cv(S_train, alpha = .8)
                     
                     list(S_train = S_train,
                          S_valid = S_valid,
                          res = OM_train)
                   },
                   .progress = TRUE)

# qs::qsave(cv_sciocv_res, "data/intermediates/230818_cv/230830_cv_sciocv_res.qs")

xx <- map(cv_scio_res,
          \(.pen) map_lgl(.pen,
                          \(.fold) isSymmetric(.fold[["res"]][["w"]])))







#~ DPM ----

cv_dpm_res <- map(unique(folds) |> set_names(),
                  \(fold){
                    S_train <- cov(mat_train[folds != fold,])
                    S_valid <- cov(mat_train[folds == fold,])
                    
                    OM_train <- DPM::dpm(S_train)
                    
                    list(S_train = S_train,
                         S_valid = S_valid,
                         res = OM_train)
                  },
                  .progress = TRUE)

# qs::qsave(cv_dpm_res, "data/intermediates/230818_cv/230829_cv_dpm_res.qs")

cv_dpm_res <- qs::qread("data/intermediates/230818_cv/230829_cv_dpm_res.qs")


xx <- map_lgl(cv_dpm_res,
              \(.fold) corpcor::make.positive.definite(.fold[["res"]]))

make.symmetric <- function(mat){
  (mat + t(mat))/2
}


cv_dpm <- cv_dpm_res |>
  enframe(name = "fold") |>
  unnest_wider(value) |>
  mutate(fold = as.integer(fold),
         res = map2(res, S_train,
                    ~ {
                      rownames(.x) <- rownames(.y)
                      colnames(.x) <- colnames(.y)
                      .x
                    }),
         res = map(res,
                   ~ .x[rowSums(is.na(.x)) < 700, colSums(is.na(.x)) < 700]),
         res = map(res, corpcor::make.positive.definite))

cv_dpm2 <- cv_dpm |>
  mutate(loss_frob = map2_dbl(S_valid, res,
                              ~ loss_frob(.x[rownames(.y),colnames(.y)],
                                          corpcor::pseudoinverse(.y))),
         loss_quad = map2_dbl(S_valid, res,
                              ~ loss_quad(.x[rownames(.y),colnames(.y)], .y)),
         loss_cv_I = map2_dbl(fold, res,
                              ~loss_tong_cv_I(mat_train[folds != .x, colnames(.y)], .y)),
         loss_cv_II = map2_dbl(fold, res,
                               ~loss_tong_cv_II(mat_train[folds != .x, colnames(.y)], .y))) |>
  mutate(adj_sparsity = map_dbl(res,
                                ~ mat_sparsity(.x)),
         adj_power_law_R2 = map_dbl(res,
                                    ~ mat_power_law(.x)))


# qs::qsave(cv_dpm2, "data/intermediates/230818_cv/230829_cv_dpm2.qs")
cv_dpm2 <- qs::qread("data/intermediates/230818_cv/230829_cv_dpm2.qs")




cv_dpm2 |>
  select(where(is.atomic))




# Source of error ----

cv_quic <- qs::qread("data/intermediates/230818_cv/230828_cv_quic.qs")


#~ Tong CV I ----


# from Tong CV I method, take y(1) as training set, y(2) as validation set
# then we want to predict validation set from training set
Y1 <- mat_train[folds == 4,]
Y2 <- mat_train[folds != 4,]

dim(Y1)
dim(Y2)

S <- cov(t(rbind(Y1,Y2)))

S11 <- S[1:nrow(Y1), 1:nrow(Y1)]
S22 <- S[(nrow(Y1) + 1):(nrow(Y1)+nrow(Y2)),(nrow(Y1) + 1):(nrow(Y1)+nrow(Y2))]
S21 <- S[(nrow(Y1) + 1):(nrow(Y1)+nrow(Y2)), 1:nrow(Y1)]

OM <- QUIC::QUIC(S, rho = .1)

OM21 <- OM$X[(nrow(Y1) + 1):(nrow(Y1)+nrow(Y2)), 1:nrow(Y1)]
OM11 <- OM$X[1:nrow(Y1), 1:nrow(Y1)]

W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
W_alternative <- solve(S22) %*% S21 # from an exact inversion (note same result with OM = solve(S))

plot(W, W_alternative)

# predict train set from validation set
est_Y1 <- t(W) %*% Y2

dim(est_Y1)

plot(Y1, est_Y1)

residuals <- est_Y1 - Y1

qqnorm(residuals)
qqline(residuals)

plot(as.numeric(residuals))
pheatmap::pheatmap(abs(residuals),
                   cluster_rows = FALSE,cluster_cols = FALSE)

xx <- colSums(abs(residuals))

head(sort(xx, decreasing = TRUE))

plot(Y1[,"B0285.1c.2"], abs(residuals[,"B0285.1c.2"]))

which(colnames(Y1) == "Y106G6H.2.1")

dim(W)

pca <- prcomp(W)
tibble(sample_id = rownames(Y2),
       PC1 = pca$x[,1],
       PC2 = pca$x[,2]) |>
  ggplot(aes(x = PC1, y = PC2)) +
  theme_classic() +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id))

pca2 <- prcomp(t(W))
tibble(sample_id = rownames(Y1),
       PC1 = pca2$x[,1],
       PC2 = pca2$x[,2]) |>
  ggplot(aes(x = PC1, y = PC2)) +
  theme_classic() +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id))


# a few really bad or many somewhat bad?
cumulated_residuals <- colSums(abs(residuals))

hist(cumulated_residuals, breaks = 150)
head(sort(cumulated_residuals, decreasing = TRUE))

pca <- prcomp(t(Y2))
tibble(id = colnames(Y2),
       PC1 = pca$x[,1],
       PC2 = pca$x[,2],
       cumul_resid = cumulated_residuals[id],
       residuals = if_else(cumul_resid > 20, "high",
                        if_else(cumul_resid < 7, "low",
                                "average"))) |>
  ggplot(aes(x = PC1, y = PC2)) +
  theme_classic() +
  geom_point(aes(color = cumul_resid, shape = residuals, size = residuals)) +
  ggrepel::geom_text_repel(aes(label = id)) +
  scale_size_manual(values = c(1,3,3))

pheatmap::pheatmap(residuals,
                   cluster_rows = FALSE, cluster_cols = FALSE)

mat_train
OM <- QUIC::QUIC(mat_train1[folds != 4,], rho = .1)
OM <- OM$X
Y <- mat_train1[folds != 4,]
loss_tong_cv_I(t(Y),OM)

loss_cv_I = map2_dbl(fold, res,
                     ~loss_tong_cv_I(mat_train[folds != .x,], .y[["X"]]))

OM <- cv_quic$res[[16]][["X"]]
dimnames(OM) <- dimnames(cv_quic$S_valid[[16]])


loss_tong_cv_I(mat_train1[folds != 4,], cv_quic$res[[16]][["X"]])


# residuals within a SF/PSI
plot(residuals[,10])
hist(abs(residuals))

#~~ re-normalizing ----

hist(mat_sf_train)

par(opar)
opar <- par(no.readonly = TRUE)

par(mfrow = c(3,3))
walk(sample(colnames(mat_sf_train2), 9),
     ~{
       hist(mat_sf_train2[,.x],
            main = .x)
     })
walk(sample(colnames(mat_psi_train2), 9),
     ~{
       hist(mat_psi_train2[,.x],
            main = .x)
     })



mat_sf_train2 <- huge::huge.npn(mat_sf_train)
mat_psi_train2 <- huge::huge.npn(mat_psi_train)


mat_train <- cbind(mat_sf_train2, mat_psi_train2)
mat_train1 <- cbind(mat_sf_train, mat_psi_train)



# transformed values
matplot(mat_train1[,1:3], mat_train[,1:3],
        # type = 'l',
        xlab = "Untransformed",
        ylab = "transformed with NPN")


sorted1 <- apply(mat_train1, 2, sort)
sorted2 <- apply(mat_train, 2, sort)

matplot(sorted1, sorted2,
        type = 'l',
        xlab = "Untransformed",
        ylab = "transformed with NPN")


#~~ rerun ----

# from Tong CV I method, take y(1) as training set, y(2) as validation set
# then we want to predict validation set from training set
Y1_npn <- mat_train[folds == 4,]
Y2_npn <- mat_train[folds != 4,]

Y1_untr <- mat_train1[folds == 4,]
Y2_untr <- mat_train1[folds != 4,]

dim(Y1_npn)
dim(Y2_npn)
dim(Y1_untr)
dim(Y2_untr)


S_npn <- cov(t(rbind(Y1_npn,Y2_npn)))
S_untr <- cov(t(rbind(Y1_untr,Y2_untr)))


OM_npn <- QUIC::QUIC(S_npn, rho = .1)
OM_untr <- QUIC::QUIC(S_untr, rho = .1)

OM21_npn <- OM_npn$X[(nrow(Y1) + 1):(nrow(Y1)+nrow(Y2)), 1:nrow(Y1)]
OM21_untr <- OM_untr$X[(nrow(Y1) + 1):(nrow(Y1)+nrow(Y2)), 1:nrow(Y1)]
OM11_npn <- OM_npn$X[1:nrow(Y1), 1:nrow(Y1)]
OM11_untr <- OM_untr$X[1:nrow(Y1), 1:nrow(Y1)]

W_npn <- - OM21_npn %*% solve(OM11_npn) # based on the estimated precision matrix
W_untr <- - OM21_untr %*% solve(OM11_untr) # based on the estimated precision matrix


# predict train set from validation set
est_Y1_npn <- t(W_npn) %*% Y2_npn
est_Y1_untr <- t(W_untr) %*% Y2_untr

dim(est_Y1_npn)
dim(est_Y1_untr)

plot(Y1_npn, est_Y1_npn)
plot(Y1_untr, est_Y1_untr)

residuals_npn <- est_Y1_npn - Y1_npn
residuals_untr <- est_Y1_untr - Y1_untr


# how are residuals higher?

plot(Y1_untr, Y1_npn)
plot(est_Y1_untr, est_Y1_npn)
plot(residuals_untr, residuals_npn)

plot(Y1_untr, residuals_untr)
plot(Y1_npn, residuals_npn)

# some examples
s <- sample(colnames(Y))


var_exp_untr <- 1 - abs(residuals_untr)/sqrt(sum(Y1_untr^2))
var_exp_npn <- 1 - abs(residuals_npn)/sqrt(sum(Y1_npn^2))

pheatmap::pheatmap(var_exp_npn,
                   cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap::pheatmap(var_exp_untr,
                   cluster_rows = FALSE, cluster_cols = FALSE)

hist(var_exp_npn, breaks = 50)
hist(var_exp_untr, breaks = 50)
#






#~ frobenius ----
loss_frob(cv_quic$S_valid[[16]], cv_quic$res[[16]][["W"]])

Sts <- cv_quic$S_valid[[16]]
Str <- cv_quic$res[[16]][["W"]]
dimnames(Str) <- dimnames(Sts)

loss_frob(Sts, Str)

dim(Str)

plot(Str, Sts)

image(abs(Sts - Str))

hist(abs(Sts - Str), breaks = 50) ; abline(v = .5, col = 'gray')
table(abs(Sts - Str) > .5)


image(abs(Sts - Str) > .5)

Sts[1:3,1:3]

nb_psi
Sts[1:3, 615:622]

pheatmap::pheatmap(1L*(abs(Sts - Str) > .5),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE)

xx <- rowSums(abs(Sts - Str) > .5)

head(sort(xx, decreasing = TRUE))

plot(Str["Y106G6H.2.1",], Sts["Y106G6H.2.1",], col = rep(c("gray", "darkred"), times = c(nb_sf, nb_psi)))

Str2 <- cov(mat_train[folds != 4,])

dat_valid <- mat_train[folds == 4,]
image(dat_valid)

plot(dat_valid[,"Y106G6H.2.1"])




#






#### Ground truth comparison ----

#~ Create consensus QUIC network ----

res_quic <- cv_quic |> filter(rho == 0.1) |> pull(res) |> lapply(\(.x) .x[["X"]][1:nb_sf, (nb_sf+1):(nb_sf+nb_psi)])

dim(res_quic[[1]])

par(opar)
opar <- par(no.readonly = TRUE)

par(mfrow = c(1,5))
walk(res_quic, image)

res_quic[[5]][1:4,1:4]

arr_quic <- array(unlist(res_quic), dim = c(nb_sf, nb_psi, length(unique(folds))))

all.equal(arr_quic[,,5], res_quic[[5]])

mat_quic <- rowMeans(arr_quic, dims = 2)

dim(mat_quic)

image(mat_quic)




#~ Compare with ground truth ----

sf_targets <- jsonlite::read_json("data/export_for_arman/sf_targets_v2.json")

sf2target <- tibble(sf_id = map_chr(sf_targets, \(x) x[["SF"]]),
                    target_id = map(sf_targets, \(x) x[["targets"]])) |>
  mutate(sf_id = if_else(sf_id == "mec-8 ad", "mec-8", sf_id),
         sf_id = map_chr(sf_id,
                         \(x) `if`(startsWith(x, "WBGene"), x, s2i(x, gids, warn_missing = TRUE)))) |>
  unnest_longer(target_id) |>
  filter(target_id != "0") |>
  mutate(sf_name = i2s(sf_id, gids, warn_missing = FALSE),
         target_name = i2s(target_id, gids, warn_missing = TRUE))


events_coordinates <- read_tsv("data/export_for_arman/221111_events_coordinates.tsv") |>
  select(event_id, target_id = gene_id)

sf_tx2g <- sf_expression |>
  select(transcript_id, gene_id, gene_name) |>
  add_row(transcript_id = "(Intercept)",
          gene_id =  "(Intercept)",
          gene_name = "(Intercept)") |>
  distinct()

convert_sf_tx2g <- function(tx_names, warn_missing = TRUE){
  res <- sf_tx2g$gene_id[match(tx_names, sf_tx2g$transcript_id, incomparables = NA)]
  if (warn_missing && any(is.na(res))) {
    warning("i2s: ", sum(is.na(res)), " tx names could not be converted. NA are returned.")
  }
  res
}


dimnames(mat_quic) <- list(colnames(mat_sf),
                           colnames(mat_psi))

exp_sf2target <- mat_quic |>
  as.data.frame() |>
  rownames_to_column("sf_name") |>
  pivot_longer(-sf_name,
               names_to = "event_id",
               values_to = "strength") |>
  filter(strength != 0) |>
  left_join(events_coordinates, by = "event_id") |>
  mutate(sf_id = convert_sf_tx2g(sf_name))


rel_exp <- exp_sf2target |> mutate(rel = paste0(sf_id,"-",target_id)) |> pull(rel) |> unique()
rel_known <- sf2target |>
  filter(target_id %in% events_coordinates$target_id) |>
  mutate(rel = paste0(sf_id,"-",target_id)) |>
  pull(rel) |> 
  unique()
rel_rand <- exp_sf2target |> mutate(target_id = sample(target_id)) |>
  mutate(rel = paste0(sf_id,"-",target_id)) |> pull(rel) |> unique()

table(rel_exp %in% rel_known)
table(rel_rand %in% rel_known)
table(rel_known %in% rel_exp)
table(rel_known %in% rel_rand)


exp_sf2target |>
  left_join(sf2target |> select(sf_id, target_id) |> filter(target_id %in% events_coordinates$target_id) |>add_column(known = 1) |> distinct(),
            by = c("sf_id", "target_id")) |>
  mutate(known = if_else(is.na(known), 0, 1) |> as.factor()) |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = known, y = abs(strength)))

exp_sf2target |>
  left_join(sf2target |> select(sf_id, target_id) |> filter(target_id %in% events_coordinates$target_id) |> add_column(known = 1) |> distinct(),
            by = c("sf_id", "target_id")) |>
  mutate(known = if_else(is.na(known), 0, 1) |> as.factor()) |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = abs(strength), color = known), adjust = 2)




### ~~~~~~~~~~ ----
r_gg <- r_gs <- r_ss <- .05
regul_mat <- rbind(cbind(matrix(r_gg, nrow=nb_sf, ncol = nb_sf),
                         matrix(r_gs, nrow=nb_sf, ncol = nb_psi)),
                   cbind(matrix(r_gs, nrow=nb_psi, ncol = nb_sf),
                         matrix(r_ss, nrow=nb_psi, ncol = nb_psi)))



quic <- QUIC::QUIC(S_train, rho = regul_mat)

# subsample for faster plotting
srows <- sample(nrow(S_train), 100)
scols <- sample(ncol(S_train), 100)

plot(S_train[srows, scols], quic$W[srows, scols],
     xlab = "covariance matrix train set",
     ylab = "Inverse of inferred precision matrix")


# Partial covariance matrices
S11 <- cov(t(Yp1))*(T-1)/T
S12 <- cov(t(Yp1), t(Yp2))*(T-1)/T
S21 <- cov(t(Yp2), t(Yp1))*(T-1)/T
S22 <- cov(t(Yp2))*(T-1)/T


# Functions to test ----


test_metrics <- function(adj, mat_test){
  nb_sf <- ncol(adj)
  nb_psi <- nrow(adj)
  
  predicted_psi <- mat_test[,1:nb_sf] %*% t(adj)
  
  mat_psi_test <- mat_test[,(nb_sf+1):(nb_sf+nb_psi)]
  nb_tests <- length(mat_psi_test)
  
  mae <- sum_na(abs(predicted_psi - mat_psi_test))/nb_tests
  sign_dpsi_matches <- sign(predicted_psi) == sign(mat_psi_test)
  accuracy <- sum_na(sign_dpsi_matches)/nb_tests
  
  list(mae = mae, accuracy = accuracy)
}

mat_metrics <- function(adj){
  
  sparsity <- sum_na(adj == 0)/length(adj)
  
  connectivity <- colSums(adj != 0)
  k <- as.numeric(names(table(connectivity)))[-1]
  p_k <- as.numeric(table(connectivity))[-1]
  R2 <-   cor(log10(k), log10(p_k))^2
  
  list(sparsity = sparsity, R2 = R2)
}



sum_na <- partial(sum, na.rm = TRUE)

f_dpm_nonreg <- function(mat_train, nb_sf, nb_psi){
  
  dpm <- DPM::dpm(mat_train)
  
  adjacency <- dpm[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]
  
  rownames(adjacency) <- colnames(mat_train)[(nb_sf+1):(nb_sf+nb_psi)]
  colnames(adjacency) <- colnames(mat_train)[1:nb_sf]
  adjacency
}


f_dpm_reg <- function(mat_train, nb_sf, nb_psi){
  
  dpm <- DPM::reg.dpm(mat_train)
  
  #  adjacency matrix
  adjacency <- dpm[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]
  
  rownames(adjacency) <- colnames(mat_train)[(nb_sf+1):(nb_sf+nb_psi)]
  colnames(adjacency) <- colnames(mat_train)[1:nb_sf]
  adjacency
}


f_arac <- function(mat_train, nb_sf, nb_psi, eps, estimator = "spearman", discretization = "none"){
  
  mim <- minet::build.mim(mat_train,
                          estimator = estimator,
                          disc = discretization)
  arac <- minet::aracne(mim, eps = eps)
  
  adjacency <- arac[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]
  adjacency
}

f_quic <- function(mat_train, nb_sf, nb_psi, r_gg = .07, r_ss = .0675, r_gs = .055){
  
  regul_mat <- rbind(cbind(matrix(r_gg, nrow=nb_sf, ncol = nb_sf),
                           matrix(r_gs, nrow=nb_sf, ncol = nb_psi)),
                     cbind(matrix(r_gs, nrow=nb_psi, ncol = nb_sf),
                           matrix(r_ss, nrow=nb_psi, ncol = nb_psi)))
  
  quic <- QUIC::QUIC(cov(mat_train), rho = regul_mat)
  
  adjacency <- quic$X[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]
  adjacency
}


# Tests ----

adj_dpm_nonreg <- f_dpm_nonreg(mat_train, nb_sf, nb_psi)
adj_dpm_reg <- f_dpm_reg(mat_train, nb_sf, nb_psi)
adj_arac <- f_arac(mat_train, nb_sf, nb_psi, eps = 0)
adj_quic <- f_quic(mat_train, nb_sf, nb_psi)


mat_metrics(adj_quic)
test_metrics(adj_quic, mat_test)


adj_arac_bin <- (adj_arac !=0)*1L
mat_metrics(adj_arac)
test_metrics(adj_arac, mat_test)
test_metrics(adj_arac, mat_train)

arac_cv_i <- function(i, eps){
  adj <- f_arac(mat_train[-i,], nb_sf, nb_psi, eps)
  test_metrics(adj, mat_train[i,, drop = FALSE])
}

arac_cv <- function(mat_train, eps, verbose = TRUE){
  if(verbose) message("With eps = ", eps)
  map(seq_len(nrow(mat_train)),
      ~ arac_cv_i(i = .x, eps),
      .progress = verbose) |>
    bind_rows()
}

# res1 <- arac_cv(mat_train, 0)
# 
# hist(res1$accuracy)
# plot(res1$mae, res1$accuracy)

# res_arac_cv_limited <- res_arac_cv

res_arac_cv <- tibble(eps = (0:10)/100,
                      res = map(eps, ~ arac_cv(mat_train, .x)))

res_arac_cv |>
  mutate(res = map(res, bind_rows)) |>
  unnest(res) |>
  pivot_longer(cols = -eps,
               names_to = "metric",
               values_to = "value") |>
  ggplot() +
  geom_jitter(aes(x = eps, y = value, color = metric))


res_arac_cv |>
  mutate(res = map(res, bind_rows)) |>
  unnest(res) |>
  summarize(mean_mae = mean(mae),
            sd_mae = sd(mae),
            mean_accuracy = mean(accuracy),
            sd_accuracy = sd(accuracy),
            .by = eps) |>
  pivot_longer(cols = -eps,
               names_sep = "_",
               names_to = c(".value", "metric"),
               values_to = "value") |>
  # filter(metric == "accuracy") |>
  ggplot(aes(x = eps, y = mean, ymin = mean-sd, ymax = mean+sd, color = metric)) +
  geom_point() +
  geom_errorbar(width = .01)


res_arac_cv |>
  mutate(res = map(res, bind_rows)) |>
  unnest(res) |>
  pivot_longer(cols = -eps,
               names_to = "metric",
               values_to = "value") |>
  filter(metric == "accuracy") |>
  ggplot() +
  geom_jitter(aes(x = eps, y = value, color = metric))




# CV ----

#https://www.stat.berkeley.edu/~bickel/BL2008-banding.pdf
# Bickel 2008: split in train/test, compute diffrence between sigma in train and test
#It may be surprising that using the sample covariance ˆ2 as the target in (24)
# works at all, since it is known to be a very noisy estimate of . It is, however, an
# unbiased estimate, and we found that even though (24) tends to overestimate the
# actual value of the risk, it gives very good results for choosing k.









#CV-I ----

dim(mat_train)

Y <- mat_train[1:4,1:4]

Y <- matrix(rnorm(60), 20,3)
Y <- apply(Y, 1, \(.y) .y-colMeans(Y)) |> t()



S <- (nrow(Y)-1)*cov(Y)/nrow(Y)

s11 <- S[1,1, drop=FALSE]
s1 <- S[-1,1, drop=FALSE]
S_ <- S[-1,-1, drop=FALSE]

all.equal(S,
          rbind(cbind(s11,t(s1)),
                cbind(s1, S_)))

denom <- as.numeric(s11 - t(s1) %*% S_ %*% s1)
om1 <- -(S_ %*% s1)/denom
om11 <- 1/denom

w1 <- -om1/om11

t(w1) %*% t(Y[,-1])
Y[,1]

t(corpcor::pseudoinverse(S_) %*% s1) %*% t(Y[,-1])

plot(Y[,1], t(w1) %*% t(Y[,-1])); abline(a=0,b=1)
plot(Y[,1], t(corpcor::pseudoinverse(S_) %*% s1) %*% Y[,-1]); abline(a=0,b=1)


OM <- DPM::reg.dpm(Y)
OM <- minet::build.mim(Y) |> minet::aracne()
# OM <- solve(S)
om11 <- OM[1,1]
om1 <- OM[-1, 1, drop=FALSE]

w1 <- -om1/om11


t(w1) %*% t(Y[,-1])
Y[,1]

plot(Y[,1], t(w1) %*% t(Y[,-1])); abline(a=0,b=1)

plot(solve(OM), S)

plot(solve(S), OM)


# CV-I bis ----

# their definition of cov
Spart <- list()
for(s in 1:nrow(Y)){
  Spart[[s]] <- t(Y[s,,drop=FALSE]) %*% Y[s,,drop=FALSE]
}

all.equal((nrow(Y)-1)*cov(Y)/nrow(Y),
          reduce(Spart, `+`)/(nrow(Y)))




Y <- matrix(rnorm(5*100, sd = 50),
            nrow = 100,
            ncol = 5)
Y <- apply(Y, 1, \(.y) .y-colMeans(Y)) |> t()

S <- (nrow(Y)-1)*cov(Y)/nrow(Y)
OM <- solve(S)

Yp1 <- Y[,1:3]
Yp2 <- Y[,4:5, drop = FALSE]


# Partial covariance matrices

Spart11 <- list()
for(s in 1:nrow(Yp1)){
  Spart11[[s]] <- t(Yp1[s,,drop=FALSE]) %*% Yp1[s,,drop=FALSE]
}
S11 <- reduce(Spart11, `+`)/(nrow(Y))

Spart22 <- list()
for(s in 1:nrow(Yp2)){
  Spart22[[s]] <- t(Yp2[s,,drop=FALSE]) %*% Yp2[s,,drop=FALSE]
}
S22 <- reduce(Spart22, `+`)/(nrow(Y))

Spart12 <- list()
for(s in 1:nrow(Yp1)){
  Spart12[[s]] <- t(Yp1[s,,drop=FALSE]) %*% Yp2[s,,drop=FALSE]
}
S12 <- reduce(Spart12, `+`)/(nrow(Y))

Spart21 <- list()
for(s in 1:nrow(Yp1)){
  Spart21[[s]] <- t(Yp2[s,,drop=FALSE]) %*% Yp1[s,,drop=FALSE]
}
S21 <- reduce(Spart21, `+`)/(nrow(Y))

rbind(
  cbind(S11, S12),
  cbind(S21, S22)
) |>
  all.equal(S)

OM21 <- -solve(S22) %*% S21 %*% solve(S11 - S12 %*% solve(S22) %*% S21)
OM11 <- solve(S11 - S12 %*% solve(S22) %*% S21)


rbind(
  cbind(OM11, t(OM21)),
  cbind(OM21, matrix(0,nrow(OM21),nrow(OM21)))
) |>
  (\(.x) abs(.x - OM) <= .Machine$double.eps)()

W <- - OM21 %*% solve(OM11)
W2 <- solve(S22) %*% S21
all.equal(W, W2)
nrow(W) == ncol(Yp2)
ncol(W) == ncol(Yp1)

t(W)
head(Yp2)

est_Yp1 <- Yp2 %*% W

head(Yp1)
head(est_Yp1)

est_Y <- list(nrow(Yp1))
for(s in 1:nrow(Yp1)){
  est_Y[[s]] <- W %*% t(Yp1[s,, drop = FALSE])
}

est_Y_all <- do.call(cbind, est_Y) |> t()
Yp2

plot(Yp1, est_Yp1)


estp1_b <- t(W) %*% t(Yp2)

measp1_b <- t(Yp1)

plot(measp1_b, estp1_b)





# CV-I with samples as rows like them ----
suppressPackageStartupMessages(library(tidyverse))

# create some data
Y_tib <- tibble(y1 = rnorm(20, sd = 50),
                y2 = 2*y1 + rnorm(20, sd = 10),
                y3 = y1*rnorm(20, sd = 0.05),
                y4 = y1/y2+y3,
                y5 = y2+y3 + rnorm(20, sd = .01)) |>
  scale(center = TRUE, scale = FALSE)

mat_train |> dim()
mat_train[1:2,1:3]
Y <- mat_train |> scale(center = TRUE, scale = FALSE)
Y |> dim()
Y[1:2,1:3]

# the rest is as above
T <- ncol(Y)
N <- nrow(Y)

# Covariance matrix
Spart <- list()
for(t in 1:T){
  Spart[[t]] <- Y[,t,drop=FALSE] %*% t(Y[,t,drop=FALSE])
}
S <- reduce(Spart, `+`)/T

all.equal(S,
          cov(t(Y))*(T-1)/T,
          check.attributes = FALSE)


N1 <- seq_len(ceiling(0.7 * nrow(mat_train)))
N2 <- setdiff(1:nrow(mat_train), N1)

Yp1 <- Y[N1,]
Yp2 <- Y[N2,]



# Partial covariance matrices
Spart11 <- list()
for(t in 1:T){
  Spart11[[t]] <- Yp1[,t,drop=FALSE] %*% t(Yp1[,t,drop=FALSE])
}
S11 <- reduce(Spart11, `+`)/T

Spart22 <- list()
for(t in 1:T){
  Spart22[[t]] <- Yp2[,t,drop=FALSE] %*% t(Yp2[,t,drop=FALSE])
}
S22 <- reduce(Spart22, `+`)/T

Spart12 <- list()
for(t in 1:T){
  Spart12[[t]] <- Yp1[,t,drop=FALSE] %*% t(Yp2[,t,drop=FALSE])
}
S12 <- reduce(Spart12, `+`)/T

Spart21 <- list()
for(t in 1:T){
  Spart21[[t]] <- Yp2[,t,drop=FALSE] %*% t(Yp1[,t,drop=FALSE])
}
S21 <- reduce(Spart21, `+`)/T


rbind(
  cbind(S11, S12),
  cbind(S21, S22)
) |>
  all.equal(S)

OM <- corpcor::pseudoinverse(S)


OM21 <- -solve(S22) %*% S21 %*% corpcor::pseudoinverse(S11 - S12 %*% solve(S22) %*% S21)
OM11 <- corpcor::pseudoinverse(S11 - S12 %*% solve(S22) %*% S21)

OM |> dim()
OM21 |> dim()
OM11 |> dim()

abs(OM[N1,N1] - OM11) |> hist()
abs(OM21 - OM[N2,N1]) |> hist()

rbind(
  cbind(OM11, t(OM21)),
  cbind(OM21, matrix(0,nrow(OM21),nrow(OM21)))
) |>
  (\(.x) abs(.x - OM) <= 1e-3)() |>
  dim()



W <- - OM21 %*% corpcor::pseudoinverse(OM11)

W_alternative <- solve(S22) %*% S21
all.equal(W, W_alternative)

W |> dim()

estimated_Yp1 <- t(W) %*% Yp2


plot(Yp1, estimated_Yp1,
     xlab = expression("True "*bold(y)^"(1)"*"(from data)"),
     ylab = expression("Estimated "*bold(y)^"(1)"))













# CV-II ----

# LDL decomposition https://en.wikipedia.org/wiki/Cholesky_decomposition
# with help from Ben Grossmann https://math.stackexchange.com/questions/4752822/equivalence-of-the-ldl-decomposition-with-an-upper-triangular-or-lower-triangula/4752910#4752910
modif_chol <- function(M){
  C <-  t(chol(K %*% M %*% K))
  stopifnot(all.equal(K %*% M %*% K, C %*% t(C)))
  
  S <- matrix(0,nrow(M),ncol(M)); diag(S) <- diag(C)
  
  L <- C %*% solve(S)
  D0 <- S^2
  
  stopifnot(all.equal(K %*% M %*% K, L %*% D0 %*% t(L)))
  
  
  # get to the same form as paper
  K <- matrix(c(0,0,1,
                0,1,0,
                1,0,0), nrow = nrow(L), byrow = TRUE)
  
  T <- t(K %*% L %*% K)
  D <- matrix(0,nrow(D0),ncol(D0)); diag(D) <- 1/rev(diag(D0))
  
  stopifnot(all.equal(M, t(T) %*% solve(D) %*% T))
  list(T = T, D = D)
}


Y <- matrix(rnorm(3*20),
            nrow = 20,
            ncol = 3)
Y <- apply(Y, 1, \(.y) .y-colMeans(Y)) |> t()

S <- cov(Y)
OM <- solve(S)

ch <- modif_chol(OM)

mat_T <- ch$T
mat_D <- ch$D

N <- ncol(Y)
err <- numeric(N)
for(n in 1:N){
  part_sum <- list(n-1)
  for(m in setdiff(1:N, n)){
    part_sum[[m]] <- mat_T[n,m] * Y[,m]
  }
  
  err[[n]] <- sum((reduce(part_sum, `+`) - Y[,n])^2)
}



plot(err, diag(D))




(Y[1,,drop = FALSE] %*% t(Y[1,,drop = FALSE]))/nrow(Y)




#



# Manual ----

mim <- minet::build.mim(mat_train)
arac <- minet::aracne(mim)

# adjacency
adj <- arac[(nb_sf + 1):(nb_sf + nb_psi), 1:nb_sf]

# test
predicted_psi <- mat_test[,1:nb_sf] %*% t(adj)

mat_psi_test <- sign(mat_test[,(nb_sf+1):(nb_sf+nb_psi)])
nb_tests <- length(mat_psi_test)


mae <- sum_na(abs(predicted_psi - mat_psi_test))/nb_tests
sign_dpsi_matches <- sign(predicted_psi) == sign(mat_psi_test)
accuracy <- sum_na(sign_dpsi_matches)/nb_tests

mae <- sum(abs(predicted_psi - mat_psi_test), na.rm = TRUE)/length(mat_psi_test)
mae

# Direct inversion ----

mat_sf_train_inv <- corpcor::pseudoinverse(mat_sf_train)
C <- mat_sf_train_inv %*% mat_psi_train

est_psi_test <- mat_sf[test_samples,] %*% C

plot(mat_psi_test,est_psi_test)


# SF clustering ----
# from "Benchmarking algorithms for gene regulatory network inference from single-cell
# transcriptomic data", Pratapa 2020
# 
# We first computed a co-expression network using the absolute value of the Pearson’s
# correlation coefficient as the weight of the edge between two genes. Next, we ran
# the Louvain clustering algorithm [6] on the predicted GRN (considering                                                                                                        all edges and their weights) and on the co-expression network to compute clusters. In order to compute clusters
#  using the Louvain algorithm, we converted the GRNs into undirected networks: if
# a method predicted both the edge (a, b) and (b, a), we used only the higher of
# the two edge weights

abscor <- abs(cor(mat_sf_train, mat_psi_train))
adj_full <- rbind(
  cbind(matrix(0, nrow = nb_sf, ncol = nb_sf), abscor),
  cbind(t(abscor), matrix(0, nrow = nb_psi, ncol = nb_psi))
)
colnames(adj_full) <- rownames(adj_full) <- c(colnames(mat_sf_train), colnames(mat_psi_train))

g <- igraph::graph_from_adjacency_matrix(adj_full, weighted = TRUE, mode = "undirected")
igraph::V(g)$type <- rep(c(TRUE, FALSE), times = c(nb_sf, nb_psi))
plot(g, layout = igraph::layout_as_bipartite)


cl <- igraph::cluster_louvain(g, resolution = 3)
igraph::sizes(cl) |> table()

plot(cl, y = g, size = 1)








# more ----

mat_psi_test <- mat_test[,(nb_sf+1):(nb_sf+nb_psi)]
# opar <- par()
par(mfrow = c(3,1), mar = c(3,0.5,0,.5))
j <- sample(rownames(mat_psi_test), 1)
plot(mat_psi_test[j,], predicted_psi[j,])
i <- sample(colnames(mat_psi_test), 1)
plot(mat_psi_test[,i], predicted_psi[,i])




# Cross-validation

nb_samples <- nrow(mat_train)
nb_folds <- 10

folds <- sample(ntile(seq_len(nb_samples), nb_folds)) |>
  as.factor()

map_dbl(levels(folds),
        \(.fold) f_arac(
          mat_train[folds != .fold,],
          nb_sf = nb_sf,
          nb_psi = nb_psi,
          mat_test = mat_train[folds == .fold,],
          eps = 0,
          estimator = "spearman",
          discretization = "none"
        )
)


map_dbl(levels(folds),
        \(.fold) f_dpm_nonreg(
          mat_train[folds != .fold,],
          nb_sf = nb_sf,
          nb_psi = nb_psi,
          mat_test = mat_train[folds == .fold,]
        )
)











