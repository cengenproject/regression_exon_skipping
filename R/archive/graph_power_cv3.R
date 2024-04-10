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




# 
# **no NPN transformation** ----
# mat_sf_train <- huge::huge.npn(mat_sf_train)
# mat_psi_train <- huge::huge.npn(mat_psi_train)
# 


#~ Assemble data ----

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








# QUIC ----


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
#      file = "data/intermediates/230918_cv_quic_nonpn.rda")


load("data/intermediates/230918_cv_quic_nonpn.rda")

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


map2_dbl(f_S_valid, f_S_train,
         loss_frob) |>
  round()









# NPN skeptic + QUIC transformation ----



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
  
  S_train <- huge::huge.npn(t(rbind(train_psi, train_sf)), npn.func = "skeptic")
  S_train <- S_train[! rowSums(is.na(S_train)) > 1, ! colSums(is.na(S_train)) > 1]
  
  res_quic <- QUIC::QUIC(S_train, rho = .1)
  OM <- res_quic[["X"]]
  dimnames(OM) <- dimnames(S_train)
  
  f_S_train[[fold]] <- res_quic[["W"]]
  
  
  f_OM[[fold]] <- OM
  
  psi_cols <- which(startsWith(colnames(OM), "SE_"))
  sf_cols <- which(! startsWith(colnames(OM), "SE_"))
  
  OM21 <- OM[sf_cols, psi_cols]
  OM11 <- OM[psi_cols, psi_cols]
  
  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
  
  sub_vald_sf <- vald_sf[rownames(vald_sf) %in% colnames(OM),]
  
  f_S_valid[[fold]] <- cov(t(rbind(vald_psi, sub_vald_sf)))
  
  # predict PSI set from SF, note transpose SF and result as we're not subsampling samples
  f_psi_estimated[[fold]] <- t(W) %*% sub_vald_sf
  
  f_psi_measured[[fold]] <- vald_psi
}

# save(f_S_train, f_S_valid, f_OM, f_psi_measured, f_psi_estimated,
#      file = "data/intermediates/230918_cv_quic_skeptic.rda")


load("data/intermediates/230918_cv_quic_nonpn.rda")

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


map2_dbl(f_S_valid, f_S_train,
         loss_frob) |>
  round()




# NPN truncated ----


# NPN transformation
mat_sf_train <- huge::huge.npn(mat_sf_train, npn.func = "truncation")
mat_psi_train <- huge::huge.npn(mat_psi_train, npn.func = "truncation")

# match rows
stopifnot(all.equal(rownames(mat_psi_train), rownames(mat_sf_train)))
mat_train <- cbind(mat_sf_train, mat_psi_train)


# finish
nb_psi <- ncol(mat_psi_train)
nb_sf <- ncol(mat_sf_train)

mat_test <- cbind(mat_sf[test_samples,], mat_psi[test_samples, ])
mat_sf_test <- mat_test[,1:nb_sf]
mat_psi_test <- mat_test[,(nb_sf+1):(nb_sf+nb_psi)]




#~ CV ----

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
#      file = "data/intermediates/230918_cv_quic_npntrunc.rda")




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


map2_dbl(f_S_valid, f_S_train,
         loss_frob) |>
  round()






# Box-Cox ----

mat_sf_train
b <- MASS::boxcox(mat_sf_train[,5]+100 ~ 1, plotit = FALSE)

# hist(mat_sf_train[,1])

lambda <- b$x[which.max(b$y)]
hist(((mat_sf_train[,5]+100)^lambda - 1)/lambda)

hist(mat_sf_train)

par(mfrow = c(3,3))
ex <- sample(colnames(mat_sf_train), 9)
walk(ex,
     ~{
       hist(sqrt(mat_sf_train[,.x]),
            main = .x)
     })



walk(sample(colnames(mat_psi_train), 9),
     ~{
       hist(mat_psi_train[,.x],
            main = .x)
     })



# mat_sf_train <- huge::huge.npn(mat_sf_train)
# mat_psi_train <- huge::huge.npn(mat_psi_train)





# QUIC untransformed vary rho ----


rho_vals <- c(10,2,1,.1,.05,.03)

f_S_train <- vector("list", length = length(rho_vals))
f_S_valid <- vector("list", length = length(rho_vals))
f_OM <- vector("list", length = length(rho_vals))
f_psi_measured <- vector("list", length = length(rho_vals))
f_psi_estimated <- vector("list", length = length(rho_vals))


fold <- 1

for(.rho in seq_along(rho_vals)){
  train_psi <- t(mat_train[folds != fold,(nb_sf+1):(nb_sf+nb_psi)])
  train_sf <- t(mat_train[folds != fold,1:nb_sf])
  
  vald_psi <- t(mat_train[folds == fold,(nb_sf+1):(nb_sf+nb_psi)])
  vald_sf <- t(mat_train[folds == fold,1:nb_sf])
  
  S_train <- cov(t(rbind(train_psi, train_sf)))
  
  
  res_quic <- QUIC::QUIC(S_train, rho = rho_vals[[.rho]])
  OM <- res_quic[["X"]]
  dimnames(OM) <- dimnames(S_train)
  
  f_S_train[[.rho]] <- res_quic[["W"]]
  
  
  f_OM[[.rho]] <- OM
  
  OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
  OM11 <- OM[1:nb_psi, 1:nb_psi]
  
  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
  
  
  f_S_valid[[.rho]] <- cov(t(rbind(vald_psi,vald_sf)))
  
  # predict PSI set from SF, note transpose SF and result as we're not subsampling samples
  f_psi_estimated[[.rho]] <- t(W) %*% vald_sf
  
  f_psi_measured[[.rho]] <- vald_psi
}


# save(f_S_train, f_S_valid, f_OM, f_psi_measured, f_psi_estimated,
#      file = "data/intermediates/230918_cv_quic_vary_rho.rda")


load("data/intermediates/230918_cv_quic_vary_rho.rda")

#~ PSI prediction ----

par(mfrow = c(2,3), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(f_psi_measured[[i]], f_psi_estimated[[i]], xlab = "",ylab = "", main = paste("rho =",rho_vals[[i]]))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(f_psi_measured[[i]])) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - f_psi_measured[[i]])

f_residuals |> lapply(abs) |> sapply(sum) |> round()


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


map2_dbl(f_S_valid, f_S_train,
         loss_frob) |>
  round()


map2_dbl(f_S_valid, f_OM,
         loss_quad) |>
  round()






# __ ----




# QUIC npn vary rho ----
#~ NPN transformation** ----
mat_sf_train <- huge::huge.npn(mat_sf_train)
mat_psi_train <- huge::huge.npn(mat_psi_train)



#~ Assemble data ----

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



#~ fit ----

rho_vals <- c(10,2,1,.1,.05,.03)

# note F-S_valid and f_psi_measured are actually all equal (since single fold)
f_S_train <- vector("list", length = length(rho_vals))
f_S_valid <- vector("list", length = length(rho_vals))
f_OM <- vector("list", length = length(rho_vals))
f_psi_measured <- vector("list", length = length(rho_vals))
f_psi_estimated <- vector("list", length = length(rho_vals))


fold <- 1

for(.rho in seq_along(rho_vals)){
  train_psi <- t(mat_train[folds != fold,(nb_sf+1):(nb_sf+nb_psi)])
  train_sf <- t(mat_train[folds != fold,1:nb_sf])
  
  vald_psi <- t(mat_train[folds == fold,(nb_sf+1):(nb_sf+nb_psi)])
  vald_sf <- t(mat_train[folds == fold,1:nb_sf])
  
  S_train <- cov(t(rbind(train_psi, train_sf)))
  
  
  res_quic <- QUIC::QUIC(S_train, rho = rho_vals[[.rho]])
  OM <- res_quic[["X"]]
  dimnames(OM) <- dimnames(S_train)
  
  f_S_train[[.rho]] <- res_quic[["W"]]
  
  
  f_OM[[.rho]] <- OM
  
  OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
  OM11 <- OM[1:nb_psi, 1:nb_psi]
  
  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
  
  
  f_S_valid[[.rho]] <- cov(t(rbind(vald_psi,vald_sf)))
  
  # predict PSI set from SF, note transpose SF and result as we're not subsampling samples
  f_psi_estimated[[.rho]] <- t(W) %*% vald_sf
  
  f_psi_measured[[.rho]] <- vald_psi
}


# save(f_S_train, f_S_valid, f_OM, f_psi_measured, f_psi_estimated,
#      file = "data/intermediates/230919_cv_quic_npn_vary_rho.rda")

load("data/intermediates/230919_cv_quic_npn_vary_rho.rda")


#~ PSI prediction ----

par(mfrow = c(2,3), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(f_psi_measured[[i]], f_psi_estimated[[i]], xlab = "",ylab = "", main = paste("rho =",rho_vals[[i]]))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(f_psi_measured[[i]])) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - f_psi_measured[[i]])

f_residuals |> lapply(abs) |> sapply(sum) |> round()


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


map2_dbl(f_S_valid, f_S_train,
         loss_frob) |>
  round()

map2_dbl(f_S_valid, f_OM,
         loss_quad) |>
  round()









# HUGE npn vary rho ----


fold <- 1

train_psi <- t(mat_train[folds != fold,(nb_sf+1):(nb_sf+nb_psi)])
train_sf <- t(mat_train[folds != fold,1:nb_sf])

vald_psi <- t(mat_train[folds == fold,(nb_sf+1):(nb_sf+nb_psi)])
vald_sf <- t(mat_train[folds == fold,1:nb_sf])

S_train <- cov(t(rbind(train_psi, train_sf)))
S_valid <- cov(t(rbind(vald_psi,vald_sf)))

res_huge <- huge::huge(S_train, method = "glasso", cov.output = TRUE)

rho_vals <- res_huge$lambda
f_OM <- res_huge$icov
f_S_train <- res_huge$cov


f_psi_estimated <- vector("list", length = length(rho_vals))

for(.rho in seq_along(res_huge$path)){
  
  OM <- f_OM[[.rho]]
  dimnames(OM) <- dimnames(S_train)
  
  OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
  OM11 <- OM[1:nb_psi, 1:nb_psi]
  
  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
  
  
  # predict PSI set from SF, note transpose SF and result as we're not subsampling samples
  f_psi_estimated[[.rho]] <- t(W) %*% vald_sf
}


# save(rho_vals, f_S_train, S_valid, f_OM, vald_psi, f_psi_estimated,
#      file = "data/intermediates/230919_cv_huge_npn_vary_rho.rda")




#~ PSI prediction ----

par(mfrow = c(2,5), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(vald_psi, f_psi_estimated[[i]], xlab = "",ylab = "",
                                   main = paste("rho =", round(rho_vals[[i]], 2)))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(vald_psi)) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - vald_psi)

f_residuals |> lapply(abs) |> sapply(sum) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- map(f_residuals,
                   frac_explained_var, vald_psi)

sapply(f_expl_var, mean) |> round(3)


f_expl_var |>
  map(enframe) |>
  imap(~add_column(.x, fold = .y)) |>
  list_rbind() |>
  summarize(val = sd(value), .by = name) |>
  arrange(desc(val))


map2_dbl(rep(list(S_valid), 10), f_S_train,
         loss_frob) |>
  round()

map2_dbl(rep(list(S_valid), 10), f_OM,
         loss_quad) |>
  round()







# __ ----


# Separate Nincl and Nexcl ----


#~ Inits ----

library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)

source("R/loss_functions.R")



quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs") |>
  filter(transcript_id != "R07E5.14.2")

#~ Prepare data ----

#~ PSI -----
mat_psi <- quantifs_filtered |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads)) |>
  select(event_id, sample_id, Nincl, Nexcl) |>
  pivot_wider(id_cols = sample_id,
              names_from = event_id,
              values_from = c(Nincl, Nexcl),
              names_vary = "slowest",
              names_glue = "{event_id}.{.value}"
  ) |>
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

#~ NPN transformation ----
mat_sf_train <- huge::huge.npn(mat_sf_train)
mat_psi_train <- huge::huge.npn(mat_psi_train)



#~ Assemble data ----

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









# QUIC ----

rho_vals <- c(10,2,1,.1,.05,.03)

# note F-S_valid and f_psi_measured are actually all equal (since single fold)
f_S_train <- vector("list", length = length(rho_vals))
f_S_valid <- vector("list", length = length(rho_vals))
f_OM <- vector("list", length = length(rho_vals))
f_psi_measured <- vector("list", length = length(rho_vals))
f_psi_estimated <- vector("list", length = length(rho_vals))


fold <- 1

for(.rho in seq_along(rho_vals)){
  train_psi <- t(mat_train[folds != fold,(nb_sf+1):(nb_sf+nb_psi)])
  train_sf <- t(mat_train[folds != fold,1:nb_sf])
  
  vald_psi <- t(mat_train[folds == fold,(nb_sf+1):(nb_sf+nb_psi)])
  vald_sf <- t(mat_train[folds == fold,1:nb_sf])
  
  S_train <- cov(t(rbind(train_psi, train_sf)))
  
  
  res_quic <- QUIC::QUIC(S_train, rho = rho_vals[[.rho]])
  OM <- res_quic[["X"]]
  dimnames(OM) <- dimnames(S_train)
  
  f_S_train[[.rho]] <- res_quic[["W"]]
  
  
  f_OM[[.rho]] <- OM
  
  OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
  OM11 <- OM[1:nb_psi, 1:nb_psi]
  
  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
  
  
  f_S_valid[[.rho]] <- cov(t(rbind(vald_psi,vald_sf)))
  
  # predict PSI set from SF, note transpose SF and result as we're not subsampling samples
  f_psi_estimated[[.rho]] <- t(W) %*% vald_sf
  
  f_psi_measured[[.rho]] <- vald_psi
}


# save(f_S_train, f_S_valid, f_OM, f_psi_measured, f_psi_estimated,
#      file = "data/intermediates/230919_cv_inclexcl_quic_npn_vary_rho.rda")

# load("data/intermediates/230919_cv_quic_npn_vary_rho.rda")


#~ PSI prediction ----

par(mfrow = c(2,3), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(f_psi_measured[[i]], f_psi_estimated[[i]], xlab = "",ylab = "", main = paste("rho =",rho_vals[[i]]))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(f_psi_measured[[i]])) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - f_psi_measured[[i]])

f_residuals |> lapply(abs) |> sapply(sum) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- map2(f_residuals, f_psi_measured,
                   frac_explained_var)

sapply(f_expl_var, mean) |> round(3)


# f_expl_var |>
#   map(enframe) |>
#   imap(~add_column(.x, fold = .y)) |>
#   list_rbind() |>
#   summarize(val = sd(value), .by = name) |>
#   arrange(desc(val))


map2_dbl(f_S_valid, f_S_train,
         loss_frob) |>
  round()

map2_dbl(f_S_valid, f_OM,
         loss_quad) |>
  round()










# HUGE npn vary rho ----


fold <- 1

train_psi <- t(mat_train[folds != fold,(nb_sf+1):(nb_sf+nb_psi)])
train_sf <- t(mat_train[folds != fold,1:nb_sf])

vald_psi <- t(mat_train[folds == fold,(nb_sf+1):(nb_sf+nb_psi)])
vald_sf <- t(mat_train[folds == fold,1:nb_sf])

S_train <- cov(t(rbind(train_psi, train_sf)))
S_valid <- cov(t(rbind(vald_psi,vald_sf)))

res_huge <- huge::huge(S_train, method = "glasso", cov.output = TRUE)

rho_vals <- res_huge$lambda
f_OM <- res_huge$icov
f_S_train <- res_huge$cov


f_psi_estimated <- vector("list", length = length(rho_vals))

for(.rho in seq_along(res_huge$path)){
  
  OM <- f_OM[[.rho]]
  dimnames(OM) <- dimnames(S_train)
  
  OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
  OM11 <- OM[1:nb_psi, 1:nb_psi]
  
  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
  
  
  # predict PSI set from SF, note transpose SF and result as we're not subsampling samples
  f_psi_estimated[[.rho]] <- t(W) %*% vald_sf
}


# save(rho_vals, f_S_train, S_valid, f_OM, vald_psi, f_psi_estimated,
#      file = "data/intermediates/230919_cv_inclexcl_huge_npn_vary_rho.rda")




#~ PSI prediction ----

par(mfrow = c(2,5), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(vald_psi, f_psi_estimated[[i]], xlab = "",ylab = "",
                                   main = paste("rho =", round(rho_vals[[i]], 2)))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(vald_psi)) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - vald_psi)

f_residuals |> lapply(abs) |> sapply(sum) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- map(f_residuals,
                  frac_explained_var, vald_psi)

sapply(f_expl_var, mean) |> round(3)


f_expl_var |>
  map(enframe) |>
  imap(~add_column(.x, fold = .y)) |>
  list_rbind() |>
  summarize(val = sd(value), .by = name) |>
  arrange(desc(val))


map2_dbl(rep(list(S_valid), 10), f_S_train,
         loss_frob) |>
  round()

map2_dbl(rep(list(S_valid), 10), f_OM,
         loss_quad) |>
  round()




# __ ----



# DPM on untransformed data ----


#~ Inits ----

library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)

source("R/loss_functions.R")



quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs") |>
  filter(transcript_id != "R07E5.14.2")

#~ Prepare data ----

#~ PSI -----
mat_psi <- quantifs_filtered |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads)) |>
  select(event_id, sample_id, Nincl, Nexcl) |>
  pivot_wider(id_cols = sample_id,
              names_from = event_id,
              values_from = c(Nincl, Nexcl),
              names_vary = "slowest",
              names_glue = "{event_id}.{.value}"
  ) |>
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

#~ NPN transformation ----
mat_sf_train <- huge::huge.npn(mat_sf_train)
mat_psi_train <- huge::huge.npn(mat_psi_train)



#~ Assemble data ----

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









# fit ----


fold <- 1

train_psi <- t(mat_train[folds != fold,(nb_sf+1):(nb_sf+nb_psi)])
train_sf <- t(mat_train[folds != fold,1:nb_sf])

vald_psi <- t(mat_train[folds == fold,(nb_sf+1):(nb_sf+nb_psi)])
vald_sf <- t(mat_train[folds == fold,1:nb_sf])


res_dpm <- DPM::reg.dpm(t(rbind(train_psi, train_sf)))


OM <- corpcor::make.positive.definite(res_dpm)
dimnames(OM) <- dimnames(S_train)

S_train <- solve(OM)



OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
OM11 <- OM[1:nb_psi, 1:nb_psi]

W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix


S_valid <- cov(t(rbind(vald_psi,vald_sf)))

# predict PSI set from SF, note transpose SF and result as we're not subsampling samples
psi_estimated <- t(W) %*% vald_sf



# save(f_S_train, f_S_valid, f_OM, f_psi_measured, f_psi_estimated,
#      file = "data/intermediates/230919_cv_inclexcl_quic_npn_vary_rho.rda")

# load("data/intermediates/230919_cv_quic_npn_vary_rho.rda")


#~ PSI prediction ----
plot(vald_psi,
     psi_estimated)

lm(as.numeric(psi_estimated) ~ as.numeric(vald_psi)) |>
  summary() |>
  (\(x) x[["adj.r.squared"]])()



# residuals
f_residuals <- list(psi_estimated - vald_psi)

f_residuals |> lapply(abs) |> sapply(sum) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- frac_explained_var(f_residuals[[1]], vald_psi)

mean(f_expl_var) |> round(3)


# f_expl_var |>
#   map(enframe) |>
#   imap(~add_column(.x, fold = .y)) |>
#   list_rbind() |>
#   summarize(val = sd(value), .by = name) |>
#   arrange(desc(val))

loss_frob(S_valid, S_train)
loss_quad(S_valid, OM)





# __ ----






# Separate npn validation sets ----

#~ Inits ----

library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)

source("R/loss_functions.R")



quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs") |>
  filter(transcript_id != "R07E5.14.2")

#~ Prepare data ----


mat_psi <- quantifs_filtered |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads)) |>
  select(event_id, sample_id, Nincl, Nexcl) |>
  pivot_wider(id_cols = sample_id,
              names_from = event_id,
              values_from = c(Nincl, Nexcl),
              names_vary = "slowest",
              names_glue = "{event_id}.{.value}"
  ) |>
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



# SF TPM
mat_sf <- sf_expression |>
  mutate(logTPM = log(TPM + 1)) |>
  select(transcript_id, sample_id, logTPM) |>
  pivot_wider(id_cols = sample_id,
              names_from = "transcript_id",
              values_from = "logTPM") |>
  column_to_rownames("sample_id") |>
  as.matrix()
mat_sf_train <- mat_sf[train_samples, ]




#~ Assemble data

# match rows
stopifnot(all.equal(rownames(mat_psi_train), rownames(mat_sf_train)))
mat_train <- cbind(mat_psi_train, mat_sf_train)


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








# QUIC ----

rho_vals <- c(10,2,1,.1,.05,.03)



fold <- 1



train_only <- mat_train[folds != fold, ] |>
  huge::huge.npn(verbose = FALSE)


S_train <- cov(train_only)


res_quic <- QUIC::QUIC(S_train, rho = 1, path = rho_vals)
f_OM <- lapply(1:dim(res_quic[["X"]])[[3]], \(i){
  OM <- res_quic[["X"]][,,i]
  dimnames(OM) <- dimnames(S_train)
  OM
  })


f_S_train <- lapply(1:dim(res_quic[["W"]])[[3]], \(i){
  S <- res_quic[["W"]][,,i]
  dimnames(S) <- dimnames(S_train)
  S
})



vald <- mat_train[folds == fold, ] |>
  huge::huge.npn(verbose = FALSE)
vald_psi <- t(vald[,1:nb_psi])
vald_sf <- t(vald[,(nb_psi+1):(nb_psi+nb_sf)])

S_valid <- cov(t(rbind(vald_psi,vald_sf)))
psi_measured <- vald_psi



f_psi_estimated <- lapply(f_OM,
       \(OM){
         OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
         OM11 <- OM[1:nb_psi, 1:nb_psi]
         
         W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
         t(W) %*% vald_sf
       })



# save(f_S_train, f_S_valid, f_OM, f_psi_measured, f_psi_estimated,
#      file = "data/intermediates/230919_cv_inclexcl_quic_sepnpn_vary_rho.rda")

# load("data/intermediates/230919_cv_inclexcl_quic_sepnpn_vary_rho.rda")


#~ PSI prediction ----

par(mfrow = c(2,3), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(vald_psi, f_psi_estimated[[i]], xlab = "",ylab = "", main = paste("rho =",rho_vals[[i]]))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

rho_vals

sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(vald_psi)) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - vald_psi)

f_residuals |> lapply(abs) |> sapply(sum) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- map(f_residuals,
                   frac_explained_var, vald_psi)

sapply(f_expl_var, mean) |> round(3)


# f_expl_var |>
#   map(enframe) |>
#   imap(~add_column(.x, fold = .y)) |>
#   list_rbind() |>
#   summarize(val = sd(value), .by = name) |>
#   arrange(desc(val))


map_dbl(f_S_train,
         ~loss_frob(S_valid, .x)) |>
  round()

map_dbl(f_OM,
        ~ loss_quad(S_valid, .x)) |>
  round()










# HUGE npn vary rho ----


fold <- 1
mat_train_train <- mat_train[folds != fold, ] |>
  huge::huge.npn(verbose = FALSE)
train_psi <- t(mat_train_train[,(nb_sf+1):(nb_sf+nb_psi)])
train_sf <- t(mat_train_train[,1:nb_sf])

mat_train_valid <- mat_train[folds == fold, ] |>
  huge::huge.npn(verbose = FALSE)
vald_psi <- t(mat_train_valid[,(nb_sf+1):(nb_sf+nb_psi)])
vald_sf <- t(mat_train_valid[,1:nb_sf])

S_train <- rbind(train_psi, train_sf) |>
  t() |>
  cov()
S_valid <- rbind(vald_psi,vald_sf) |>
  t() |>
  cov()

res_huge <- huge::huge(S_train, method = "glasso", cov.output = TRUE)

rho_vals <- res_huge$lambda
f_OM <- res_huge$icov
f_S_train <- res_huge$cov


f_psi_estimated <- lapply(f_OM,
                          \(OM){
                            OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
                            OM11 <- OM[1:nb_psi, 1:nb_psi]
                            
                            W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                            t(W) %*% vald_sf
                          })



# save(rho_vals, f_S_train, S_valid, f_OM, vald_psi, f_psi_estimated,
#      file = "data/intermediates/230919_cv_inclexcl_huge_npn_vary_rho.rda")




#~ PSI prediction ----

par(mfrow = c(2,5), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(vald_psi, f_psi_estimated[[i]], xlab = "",ylab = "",
                                   main = paste("rho =", round(rho_vals[[i]], 2)))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

rho_vals |> round(2)

sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(vald_psi)) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - vald_psi)

f_residuals |> lapply(abs) |> sapply(sum, na.rm = TRUE) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- map(f_residuals,
                  frac_explained_var, vald_psi)

sapply(f_expl_var, mean) |> round(3)


# f_expl_var |>
#   map(enframe) |>
#   imap(~add_column(.x, fold = .y)) |>
#   list_rbind() |>
#   summarize(val = sd(value), .by = name) |>
#   arrange(desc(val))


map2_dbl(rep(list(S_valid), 10), f_S_train,
         loss_frob) |>
  round()

map2_dbl(rep(list(S_valid), 10), f_OM,
         loss_quad) |>
  round()







# __ ----






# 4-way separate npn validation sets ----

#~ Inits ----

library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)

source("R/loss_functions.R")



quantifs_filtered <- qs::qread("data/intermediates/simultation/230206_preprocessed_quantifs_filtered.qs")
sf_expression <- qs::qread("data/intermediates/simultation/230206_preprocessed_sf_expression.qs") |>
  filter(transcript_id != "R07E5.14.2")

#~ Prepare data ----

#~ PSI -----
mat_psi <- quantifs_filtered |>
  mutate(Nincl = round(PSI * nb_reads),
         Nexcl = round((1-PSI) * nb_reads)) |>
  select(event_id, sample_id, Nincl, Nexcl) |>
  pivot_wider(id_cols = sample_id,
              names_from = event_id,
              values_from = c(Nincl, Nexcl),
              names_vary = "slowest",
              names_glue = "{event_id}.{.value}"
  ) |>
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

# *no imputation * ----


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




#~ Assemble data ----

# match rows
stopifnot(all.equal(rownames(mat_psi_train), rownames(mat_sf_train)))
mat_train <- cbind(mat_psi_train, mat_sf_train)


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








# QUIC ----

rho_vals <- c(10,2,1,.1,.05,.03)



fold <- 1



train_only <- mat_train[folds != fold, ]
train_psi <- train_only[,1:nb_psi] |>
  huge::huge.npn(verbose = FALSE) |>
  t()
train_sf <- train_only[,(nb_psi+1):(nb_psi+nb_sf)] |>
  huge::huge.npn(verbose = FALSE) |>
  t()



S_train <- cov(t(rbind(train_psi,train_sf)))


res_quic <- QUIC::QUIC(S_train, rho = 1, path = rho_vals)
f_OM <- lapply(1:dim(res_quic[["X"]])[[3]], \(i){
  OM <- res_quic[["X"]][,,i]
  dimnames(OM) <- dimnames(S_train)
  OM
})


f_S_train <- lapply(1:dim(res_quic[["W"]])[[3]], \(i){
  S <- res_quic[["W"]][,,i]
  dimnames(S) <- dimnames(S_train)
  S
})



vald <- mat_train[folds == fold, ]
vald_psi <- vald[,1:nb_psi] |>
  huge::huge.npn(verbose = FALSE) |>
  t()
vald_sf <- vald[,(nb_psi+1):(nb_psi+nb_sf)] |>
  huge::huge.npn(verbose = FALSE) |>
  t()

S_valid <- cov(t(rbind(vald_psi,vald_sf)))


f_psi_estimated <- lapply(f_OM,
                          \(OM){
                            OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
                            OM11 <- OM[1:nb_psi, 1:nb_psi]
                            
                            W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                            t(W) %*% vald_sf
                          })



# save(f_S_train, f_S_valid, f_OM, f_psi_measured, f_psi_estimated,
#      file = "data/intermediates/230919_cv_inclexcl_quic_sepnpn_vary_rho.rda")

# load("data/intermediates/230919_cv_inclexcl_quic_sepnpn_vary_rho.rda")


#~ PSI prediction ----

par(mfrow = c(2,3), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(vald_psi, f_psi_estimated[[i]], xlab = "",ylab = "", main = paste("rho =",rho_vals[[i]]))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

rho_vals

sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(vald_psi)) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - vald_psi)

f_residuals |> lapply(abs) |> sapply(sum) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- map(f_residuals,
                  frac_explained_var, vald_psi)

sapply(f_expl_var, mean) |> round(3)


# f_expl_var |>
#   map(enframe) |>
#   imap(~add_column(.x, fold = .y)) |>
#   list_rbind() |>
#   summarize(val = sd(value), .by = name) |>
#   arrange(desc(val))


map_dbl(f_S_train,
        ~loss_frob(S_valid, .x)) |>
  round()

map_dbl(f_OM,
        ~ loss_quad(S_valid, .x)) |>
  round()










# HUGE npn vary rho ----


fold <- 1



train_only <- mat_train[folds != fold, ]
train_psi <- train_only[,1:nb_psi] |>
  huge::huge.npn(verbose = FALSE) |>
  t()
train_sf <- train_only[,(nb_psi+1):(nb_psi+nb_sf)] |>
  huge::huge.npn(verbose = FALSE) |>
  t()



S_train <- cov(t(rbind(train_psi,train_sf)))





res_huge <- huge::huge(S_train, method = "glasso", cov.output = TRUE)

rho_vals <- res_huge$lambda
f_OM <- res_huge$icov
f_S_train <- res_huge$cov



vald <- mat_train[folds == fold, ]
vald_psi <- vald[,1:nb_psi] |>
  huge::huge.npn(verbose = FALSE) |>
  t()
vald_sf <- vald[,(nb_psi+1):(nb_psi+nb_sf)] |>
  huge::huge.npn(verbose = FALSE) |>
  t()
S_valid <- cov(t(rbind(vald_psi,vald_sf)))


f_psi_estimated <- lapply(f_OM,
                          \(OM){
                            OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
                            OM11 <- OM[1:nb_psi, 1:nb_psi]
                            
                            W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                            t(W) %*% vald_sf
                          })



# save(rho_vals, f_S_train, S_valid, f_OM, vald_psi, f_psi_estimated,
#      file = "data/intermediates/230919_cv_inclexcl_huge_npn_vary_rho.rda")




#~ PSI prediction ----

par(mfrow = c(2,5), mar = c(2,2,2,2))
for(i in seq_along(rho_vals)) plot(vald_psi, f_psi_estimated[[i]], xlab = "",ylab = "",
                                   main = paste("rho =", round(rho_vals[[i]], 2)))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

rho_vals |> round(2)

sapply(seq_along(rho_vals),
       \(i){
         lm(as.numeric(f_psi_estimated[[i]]) ~ as.numeric(vald_psi)) |>
           summary() |>
           (\(x) x[["adj.r.squared"]])()
       }
) |> round(3)


# residuals
f_residuals <- lapply(seq_along(rho_vals),
                      \(i) f_psi_estimated[[i]] - vald_psi)

f_residuals |> lapply(abs) |> sapply(sum, na.rm = TRUE) |> round()


# Fraction explained variance
# for each row (=event), we take the SSresidual and the variance (SStotal)
f_expl_var <- map(f_residuals,
                  frac_explained_var, vald_psi)

sapply(f_expl_var, mean) |> round(3)


# f_expl_var |>
#   map(enframe) |>
#   imap(~add_column(.x, fold = .y)) |>
#   list_rbind() |>
#   summarize(val = sd(value), .by = name) |>
#   arrange(desc(val))


map2_dbl(rep(list(S_valid), 10), f_S_train,
         loss_frob) |>
  round()

map2_dbl(rep(list(S_valid), 10), f_OM,
         loss_quad) |>
  round()




