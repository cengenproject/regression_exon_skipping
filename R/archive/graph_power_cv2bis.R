# New iteration on the graph power approach (2023-09-14)



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


# Train/test split
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




# Assemble data ----

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









#~~~~~~~~~~~~ ----



# QUIC entire matrix ----
mat_train_npn <- huge::huge.npn(mat_train)
res_quic <- QUIC::QUIC(cov(mat_train_npn), rho = .1)

OM <- res_quic$X
Strain <- res_quic$W






# Characterize precision matrix ----

#~ By Tong CV I ----


# from Tong CV I method, we try yo predict PSI from SF:
# take y(1) as PSI set, y(2) as SF
# Note: to stay consistent with paper recomputing OM with PSI in front of SF
Y_psi <- t(mat_train_npn[,(nb_sf+1):(nb_sf+nb_psi)])
Y_sf <- t(mat_train_npn[,1:nb_sf])


dim(Y_psi)
dim(Y_sf)


S <- cov(t(rbind(Y_psi,Y_sf)))
OM <- QUIC::QUIC(S, rho = .1)[["X"]]
dimnames(OM) <- dimnames(S)

OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
OM11 <- OM[1:nb_psi, 1:nb_psi]

W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix

# predict PSI set from SF, note transpose SF and result as we're not subsampling samples
est_Y_psi <- t(W) %*% Y_sf

dim(est_Y_psi)


# pca
dim(W)
pheatmap::pheatmap(W,
                   cluster_rows = FALSE,cluster_cols = FALSE)


pca <- prcomp(W)
tibble(sf_id = rownames(Y_sf),
       PC1 = pca$x[,1],
       PC2 = pca$x[,2]) |>
  ggplot(aes(x = PC1, y = PC2)) +
  theme_classic() +
  # ggrepel::geom_text_repel(aes(label = sf_id)) +
  geom_point()

pca2 <- prcomp(t(W))
tibble(event_id = rownames(Y_psi),
       PC1 = pca2$x[,1],
       PC2 = pca2$x[,2]) |>
  ggplot(aes(x = PC1, y = PC2)) +
  theme_classic() +
  # ggrepel::geom_text_repel(aes(label = event_id)) +
  geom_point()

# checking prediction quality
plot(Y_psi, est_Y_psi,
     xlab = "PSI measured",
     ylab = "Estimated PSI (from precision matrix)")

lm(as.numeric(est_Y_psi) ~ as.numeric(Y_psi)) |> summary()

residuals <- est_Y_psi - Y_psi

qqnorm(residuals)
qqline(residuals)

plot(Y_psi, residuals, xlab = "PSI measured")

plot(as.numeric(residuals))
hist(residuals, breaks = 50)

pheatmap::pheatmap(abs(residuals),
                   cluster_rows = FALSE,cluster_cols = FALSE)



# residuals within a SF/PSI
plot(Y_psi[10,], residuals[10,])


# Explained variance
expl_var <- frac_explained_var(residuals, Y_psi)
plot(rowMeans(Y_psi), expl_var, xlab = "Mean PSI measured", ylab = "Explained variance")

plot(matrixStats::rowVars(Y_psi), expl_var, xlab = "Variance PSI measured", ylab = "Explained variance")

mean(expl_var)







# In a CV setting ----

# we use one fold as validation, and the rest of folds as training.
# we use the training to compute OM and W, and look for error on validation


#~ Compute OM and W from training data ----

# from Tong CV I method, we try yo predict PSI from SF:
# take y(1) as PSI set, y(2) as SF
# Note: to stay consistent with paper recomputing OM with PSI in front of SF
train_only <- mat_train[folds != 1,] |>
  huge::huge.npn()
Y_psi <- t(train_only[,(nb_sf+1):(nb_sf+nb_psi)])
Y_sf <- t(train_only[,1:nb_sf])

dim(Y_psi)
dim(Y_sf)

S_train <- cov(t(rbind(Y_psi,Y_sf)))
res_quic <- QUIC::QUIC(S_train, rho = .1)
OM <- res_quic[["X"]]
dimnames(OM) <- dimnames(S_train)

OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
OM11 <- OM[1:nb_psi, 1:nb_psi]

W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix


#~ check on validation data ----
vald_only <- mat_train[folds == 1, ] |>
  huge::huge.npn()
vald_psi <- t(vald_only[, (nb_sf+1):(nb_sf+nb_psi)])
vald_sf <- t(vald_only[, 1:nb_sf])

dim(vald_psi)
dim(vald_sf)

# predict PSI set from SF, note transpose SF and result as we're not subsampling samples
est_vald_psi <- t(W) %*% vald_sf

dim(est_vald_psi)


# pca
dim(W)
pheatmap::pheatmap(W,
                   cluster_rows = FALSE,cluster_cols = FALSE)


pca <- prcomp(W)
tibble(sf_id = rownames(Y_sf),
       PC1 = pca$x[,1],
       PC2 = pca$x[,2]) |>
  ggplot(aes(x = PC1, y = PC2)) +
  theme_classic() +
  # ggrepel::geom_text_repel(aes(label = sf_id)) +
  geom_point()

pca2 <- prcomp(t(W))
tibble(event_id = rownames(Y_psi),
       PC1 = pca2$x[,1],
       PC2 = pca2$x[,2]) |>
  ggplot(aes(x = PC1, y = PC2)) +
  theme_classic() +
  # ggrepel::geom_text_repel(aes(label = event_id)) +
  geom_point()


# checking prediction quality
plot(vald_psi, est_vald_psi,
     xlab = "PSI measured (validation)",
     ylab = "Estimated PSI (from precision matrix on training)")
lm(as.numeric(est_vald_psi) ~ as.numeric(vald_psi)) |> summary()

residuals <- est_vald_psi - vald_psi

qqnorm(residuals)
qqline(residuals)

plot(vald_psi, residuals, xlab = "PSI measured")

plot(as.numeric(residuals))
hist(residuals, breaks = 50)

pheatmap::pheatmap(abs(residuals),
                   cluster_rows = FALSE,cluster_cols = FALSE)


# residuals within a SF/PSI
plot(vald_psi[10,], residuals[10,])
plot(vald_psi[10,], est_vald_psi[10,])


# Fraction Explained Variance

expl_var <- frac_explained_var(residuals, vald_psi)
plot(rowMeans(vald_psi), expl_var, xlab = "Mean PSI measured", ylab = "Explained variance")

plot(matrixStats::rowVars(vald_psi), expl_var, xlab = "Variance PSI measured", ylab = "Explained variance")
mean(expl_var)



#~ Compare to pseudoinverse ----

OM_pseudo <- corpcor::pseudoinverse(S_train)
dimnames(OM_pseudo) <- dimnames(S_train)

OM21 <- OM_pseudo[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
OM11 <- OM_pseudo[1:nb_psi, 1:nb_psi]

W_pseudo <- - OM21 %*% corpcor::pseudoinverse(OM11) # based on the estimated precision matrix


# predict validation PSI set from SF, note transpose SF and result as we're not subsampling samples
est_vald_psi_pseudo <- t(W_pseudo) %*% vald_sf



plot(vald_psi, est_vald_psi,
     xlab = "PSI measured (validation)",
     ylab = "Estimated PSI (from precision matrix on training)")
plot(vald_psi, est_vald_psi_pseudo,
     xlab = "PSI measured (validation)",
     ylab = "Estimated PSI (from pseudoinv matrix on training)")

lm(as.numeric(est_vald_psi_pseudo) ~ as.numeric(vald_psi)) |> summary()


plot(residuals, est_vald_psi_pseudo - vald_psi,
     xlab = "Residuals (sparse precision matrix)",
     ylab = "Residuals (pseudoinverse)")


resid_pseudo <- est_vald_psi_pseudo - vald_psi

qqnorm(resid_pseudo)
qqline(resid_pseudo)

# Fraction Explained Variance

expl_var_pseudo <- frac_explained_var(resid_pseudo, vald_psi)
plot(rowMeans(vald_psi), expl_var_pseudo, xlab = "Mean PSI measured", ylab = "Explained variance")

plot(matrixStats::rowVars(vald_psi), expl_var, xlab = "Variance PSI measured", ylab = "Explained variance")
mean(expl_var)


#~ Frobenius compare cov matrices ----


S_valid <- cov(t(rbind(vald_psi,vald_sf)))

plot(S_valid, S_train)

raw_cov_resid <- S_train - S_valid

S_hat_train <- res_quic$W

plot(S_valid, S_hat_train)

plot(S_valid[1:nb_psi,1:nb_psi], S_hat_train[1:nb_psi,1:nb_psi],
     main = "PSI-PSI", xlab = expression(S[validation]), ylab = expression(hat(S)[train]))
plot(S_valid[(nb_psi+1):(nb_psi+nb_sf),(nb_psi+1):(nb_psi+nb_sf)], S_hat_train[(nb_psi+1):(nb_psi+nb_sf),(nb_psi+1):(nb_psi+nb_sf)],
     main = "SF-SF", xlab = expression(S[validation]), ylab = expression(hat(S)[train]))
plot(S_valid[(nb_psi+1):(nb_psi+nb_sf),1:nb_psi], S_hat_train[(nb_psi+1):(nb_psi+nb_sf),1:nb_psi],
     main = "SF-PSI", xlab = expression(S[validation]), ylab = expression(hat(S)[train]))

plot(S_valid[1:nb_psi,(nb_psi+1):(nb_psi+nb_sf)], S_hat_train[1:nb_psi,(nb_psi+1):(nb_psi+nb_sf)],
     main = "PSI-SF", xlab = expression(S[validation]), ylab = expression(hat(S)[train]))

plot(diag(S_valid), diag(S_hat_train))

xx <- which(S_hat_train > 1.4)
length(xx)

pheatmap::pheatmap(1*(S_hat_train > 1.4),
                   cluster_rows = FALSE, cluster_cols = FALSE)

cov_resid <- S_hat_train - S_valid


hist(raw_cov_resid, breaks = 50)
hist(cov_resid, breaks = 50, add = TRUE, col = rgb(.8,0,.2,.5))

pheatmap::pheatmap(cov_resid,
                   cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap::pheatmap(raw_cov_resid,
                   cluster_rows = FALSE, cluster_cols = FALSE)

loss_frob(Sts = S_valid, Str = S_hat_train)



#~ Quadratic compare OMtrain to Svalid ----

prod <- S_valid %*% OM

pheatmap::pheatmap(prod - diag(1, nrow = nrow(OM)),
                   cluster_rows = FALSE, cluster_cols = FALSE)

hist(prod - diag(1,nrow(OM)), breaks = 50)

hist(diag(prod))

prod_raw <- S_valid %*% corpcor::pseudoinverse(S_train)
colnames(prod_raw) <- rownames(prod_raw)
hist(diag(prod_raw))

pheatmap::pheatmap(prod_raw,
                   cluster_rows = FALSE, cluster_cols = FALSE)

hist(prod - diag(1,nrow(OM)), breaks = 50)

plot(prod_raw, prod)



#~ Permutation test ----


Y_psi <- t(mat_train[folds != 1,(nb_sf+1):(nb_sf+nb_psi)])
Y_sf <- t(mat_train[folds != 1,1:nb_sf])

permut_est_vald_psi <- replicate(20,
                                 {
                                   Y_sf_perm <- t(apply(Y_sf, 1, sample))
                                   S_train <- cov(t(rbind(Y_psi,Y_sf_perm)))
                                   res_quic <- QUIC::QUIC(S_train, rho = .1)
                                   OM <- res_quic[["X"]]
                                   
                                   OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
                                   OM11 <- OM[1:nb_psi, 1:nb_psi]
                                   
                                   W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                                   
                                   vald_sf <- t(mat_train[folds == 1, 1:nb_sf])
                                   
                                   t(W) %*% vald_sf
                                 }
)



par(mfrow = c(1,3))
walk(1:3, ~plot(vald_psi, permut_est_vald_psi[,,.x]))
par(mfrow = c(1,1))

map_dbl(1:dim(permut_est_vald_psi)[[3]],
        ~ summary(lm(as.numeric(permut_est_vald_psi[,,.x]) ~ as.numeric(vald_psi)))[["adj.r.squared"]])


perm_resid <- map(1:dim(permut_est_vald_psi)[[3]],
                  ~ permut_est_vald_psi[,,.x] - vald_psi)

hist(unlist(perm_resid))

pheatmap::pheatmap(perm_resid[[1]], cluster_rows = FALSE, cluster_cols = FALSE)

plot(residuals, perm_resid[[1]],
     xlab = "Residuals (non-permuted)", ylab = "Residuals (permuted)")

plot(vald_psi, perm_resid[[7]],
     xlab = "Measured PSI", ylab = "Residuals (permuted)")


perm_sum_res <- map_dbl(perm_resid,
                        ~ sum(abs(.x)))
obs_sum_res <- sum(abs(residuals))

mean(obs_sum_res >= perm_sum_res)

hist(perm_sum_res, xlim = c(1000, 3000),
     xlab = "Sum of absolute residuals", main = NULL)
abline(v = obs_sum_res, col = 'darkred')





#







# Cross validation ----

# to estimate what part of error comes from particular split, use more folds

fold_names <- sort(unique(folds))

f_S_train <- vector("list", length = length(fold_names))
f_S_valid <- vector("list", length = length(fold_names))
f_OM <- vector("list", length = length(fold_names))
f_psi_measured <- vector("list", length = length(fold_names))
f_psi_estimated <- vector("list", length = length(fold_names))


for(fold in fold_names){
  train_psi <- t(mat_train[folds != fold,(nb_sf+1):(nb_sf+nb_psi)])
  train_sf <- t(mat_train[folds != fold,1:nb_sf])
  
  vald_psi <- t(mat_train[folds == fold,(nb_sf+1):(nb_sf+nb_psi)])
  vald_sf <- t(mat_train[folds == fold,1:nb_sf])
  
  S_train <- cov(t(rbind(train_psi, train_sf)))
  
  
  res_quic <- QUIC::QUIC(S_train, rho = .1)
  OM <- res_quic[["X"]]
  dimnames(OM) <- dimnames(S_train)
  
  f_S_train[[fold]] <- res_quic[["W"]]
  
  
  f_OM[[fold]] <- OM
  
  OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
  OM11 <- OM[1:nb_psi, 1:nb_psi]
  
  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
  
  
  f_S_valid[[fold]] <- cov(t(rbind(vald_psi,vald_sf)))
  
  # predict PSI set from SF, note transpose SF and result as we're not subsampling samples
  f_psi_estimated[[fold]] <- t(W) %*% vald_sf
  
  f_psi_measured[[fold]] <- vald_psi
}

# save(f_S_train, f_S_valid, f_OM, f_psi_measured, f_psi_estimated,
#      file = "data/intermediates/230918_cv_quic.rda")


load("data/intermediates/230918_cv_quic.rda")

#~ PSI prediction ----

par(mfrow = c(2,3), mar = c(2,2,2,2))
for(i in fold_names) plot(f_psi_measured[[i]], f_psi_estimated[[i]], xlab = "",ylab = "")
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


sapply(fold_names,
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(f_psi_measured[[i]])) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(fold_names,
                      \(i) f_psi_estimated[[i]] - f_psi_measured[[i]])

f_residuals |> sapply(abs) |> sapply(sum) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- map2(f_residuals, f_psi_measured,
                   frac_explained_var)

sapply(f_expl_var, mean) |> round(3)


f_expl_var |>
  map(enframe) |>
  imap(~add_column(.x, fold = .y)) |>
  list_rbind() |>
  summarize(val = sd(value), .by = name) |>
  arrange(desc(val))






#~ Frobenius compare cov matrices ----



plot(f_S_valid[[1]], f_S_train[[1]])

f_cov_resid <- map2(f_S_train, f_S_valid,
                    ~ .y - .x)

f_cov_resid |>
  map(as.numeric) |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  list_rbind() |>
  mutate(fold = as.factor(fold)) |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = value, fill = fold), alpha = .2) +
  xlab("Residual (covariance estimated from training vs measured in validation)")




hist(unlist(f_cov_resid), breaks = 50)


pheatmap::pheatmap(f_cov_resid[[5]],
                   cluster_rows = FALSE, cluster_cols = FALSE)

map2_dbl(f_S_valid, f_S_train,
         loss_frob) |>
  round()



#~ Quadratic compare OMtrain to Svalid ----

f_prod <- map2(f_S_valid, f_OM,
               ~ .x %*% .y)


for(i in seq_along(f_prod)) pheatmap::pheatmap(f_prod[[i]] - diag(1, nrow = nrow(f_OM[[i]])),
                                               cluster_rows = FALSE, cluster_cols = FALSE,
                                               filename = paste0("data/intermediates/230918_heatmap_quad/", i, ".png"))

f_prod |>
  map(~ .x - diag(1, nrow = nrow(f_OM[[1]]))) |>
  map(as.numeric) |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  list_rbind() |>
  mutate(fold = as.factor(fold)) |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = value, fill = fold), alpha = .2) +
  xlab(expression(paste("Residual (",S[validation], "*", hat(Omega)[train] - Id, ")")))


f_quad_resid <- map(f_prod, ~ .x - diag(1, nrow = nrow(f_OM[[1]])))

head(colSums(abs(f_quad_resid[[1]])))
head(rowSums(abs(f_quad_resid[[1]])))

hist(colSums(abs(f_quad_resid[[1]])))

sum_abs_resid <- f_prod |>
  map(~ .x - diag(1, nrow = nrow(f_OM[[1]]))) |>
  map(abs) |>
  map(colSums)

sum_abs_resid |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  list_rbind() |>
  mutate(fold = as.factor(fold)) |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = value, fill = fold), alpha = .2) +
  xlab("Sum of absolute residuals per column")


sum_abs_resid |>
  map_dbl(max) |>
  round()

# is it always the same column, how stable are these colSum residuals?
sum_abs_resid |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  list_rbind() |>
  mutate(fold = as.factor(fold)) |>
  arrange(desc(value)) |>
  mutate(event_id = fct_inorder(name)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = event_id, y = value, color = fold))


plot(sum_abs_resid[[4]],
     sum_abs_resid[[5]],
     xlab = "Sum of absolute residuals in fold 4",
     ylab = "Sum of absolute residuals in fold 5")

plot(sum_abs_resid[[2]], sum_abs_resid[[3]],
     xlab = "Sum of absolute residuals in fold 2",
     ylab = "Sum of absolute residuals in fold 3")



sum_abs_resid |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  map(~ arrange(.x, desc(value))) |>
  map(slice_head, n = 5) |>
  map(pull, name) |>
  set_names(paste0("fold_", 1:5)) |>
  UpSetR::fromList() |>
  UpSetR::upset()


sum_abs_resid |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  map(~ arrange(.x, desc(value))) |>
  map(slice_head, n = 5)
f_prod |>
  map(~ .x - diag(1, nrow = nrow(f_OM[[1]]))) |>
  map(abs) |>
  map(colSums) |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  map(~ arrange(.x, desc(value))) |>
  map(slice_tail, n = 5)

all_res_means <- sum_abs_resid |>
  map_dbl(mean)

all_res_sds <- sum_abs_resid |>
  map_dbl(sd)
hist(sum_abs_resid[[1]], breaks = 50)
abline(v = all_res_means[[1]] - 2*all_res_sds[[1]], col = 'darkred')
abline(v = all_res_means[[1]] + 2*all_res_sds[[1]], col = 'darkred')

extreme_columns <- vector("list", length(sum_abs_resid))
for(i in seq_along(sum_abs_resid)){
  extr <- which(sum_abs_resid[[i]] < all_res_means[[i]] - 2*all_res_sds[[i]] | 
                  sum_abs_resid[[i]] > all_res_means[[i]] + 2*all_res_sds[[i]])
  extreme_columns[[i]] <- names(sum_abs_resid[[i]])[extr]
}

extreme_columns |>
  setNames(paste0("fold_", 1:5)) |>
  UpSetR::fromList() |>
  UpSetR::upset()

most_extreme <- reduce(extreme_columns, intersect)

sum_abs_resid |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  list_rbind() |>
  filter(name %in% most_extreme) |>
  pivot_wider(id_cols = name,
              values_from = "value",
              names_from = "fold")

not_most_extreme <- names(sum_abs_resid[[1]]) |>
  setdiff(unlist(extreme_columns)) |>
  sample(5)

sum_abs_resid |>
  map(enframe) |>
  imap(~ add_column(.x, fold = .y)) |>
  list_rbind() |>
  filter(name %in% not_most_extreme) |>
  pivot_wider(id_cols = name,
              values_from = "value",
              names_from = "fold")





























