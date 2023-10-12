# Based on graph_power_cv and graph_power_cv3; only NPN/imputation within each fold, test methods


# Inits ----

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






# ********** ----







# QUIC ----

rho_vals <- c(10, 5, 2, 1, .5, .1, .05, .03) |> set_names()
rho_vals <- c(.4, .3, .2, .09, .08, .07, .06) |> set_names()

fold_names <- sort(unique(folds)) |> set_names()



res_quic <- map(fold_names,
                \(fold){
                  
                  train_only <- mat_train[folds != fold, ] |>
                    huge::huge.npn(verbose = FALSE)
                  train_psi <- t(train_only[,1:nb_psi])
                  train_sf <- t(train_only[,(nb_psi+1):(nb_psi+nb_sf)])
                  
                  S_train <- cov(t(rbind(train_psi,train_sf)))
                  
                  # run estimation!
                  res_quic <- QUIC::QUIC(S_train, rho = 1, path = rho_vals)
                  
                  
                  # Extract results
                  f_OM <- lapply(1:dim(res_quic[["X"]])[[3]], \(i){
                    OM <- res_quic[["X"]][,,i]
                    dimnames(OM) <- dimnames(S_train)
                    OM
                  }) |> set_names(rho_vals)
                  
                  f_S_train <- lapply(1:dim(res_quic[["W"]])[[3]], \(i){
                    S <- res_quic[["W"]][,,i]
                    dimnames(S) <- dimnames(S_train)
                    S
                  }) |> set_names(rho_vals)
                  
                  # get validation set
                  vald_only <- mat_train[folds == fold, ] |>
                    huge::huge.npn(verbose = FALSE)
                  vald_psi <- t(vald_only[,1:nb_psi])
                  vald_sf <- t(vald_only[,(nb_psi+1):(nb_psi+nb_sf)])
                  
                  S_valid <- cov(t(rbind(vald_psi,vald_sf)))
                  
                  
                  f_psi_estimated <- lapply(f_OM,
                                            \(OM){
                                              OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
                                              OM11 <- OM[1:nb_psi, 1:nb_psi]
                                              
                                              W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                                              t(W) %*% vald_sf
                                            }) |> set_names(rho_vals)
                  
                  
                  list(S_train = f_S_train,
                       OM = f_OM,
                       S_valid = S_valid,
                       psi_measured = vald_psi,
                       psi_estimated = f_psi_estimated)
                },
                .progress = TRUE)

# qs::qsave(res_quic, "data/intermediates/230920_cv/230921_quic.qs")

# previously ran in 2 passes, combined below
# res_quic <- qs::qread("data/intermediates/230920_cv/230920_quic.qs")
# res_quic <- qs::qread("data/intermediates/230920_cv/230921_quic.qs")

# add names a posteriori, only for first pass
res_quic2 <- map(res_quic,
                 ~ {
                   .x$OM <- .x$OM |> set_names(rho_vals)
                   .x$S_train <- .x$S_train |> set_names(rho_vals)
                   .x$psi_estimated <- .x$psi_estimated |> set_names(rho_vals)
                   .x
                 })


#~ PSI prediction ----

par(mfrow = c(2,4), mar = c(2,2,2,2))
for(i in seq_along(res_quic2$`1`$OM)) plot(res_quic2$`1`$psi_measured, res_quic2$`1`$psi_estimated[[i]], xlab = "",ylab = "", main = paste("rho =",rho_vals[[i]]))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


# extract values
res_quict <- transpose(res_quic2)

tib_quic <- expand_grid(penalty = rho_vals, fold = fold_names) |>
  mutate(S_train = map2(fold, penalty,
                        ~ res_quict$S_train[[as.character(.x)]][[as.character(.y)]]),
         OM = map2(fold, penalty,
                   ~ res_quict$OM[[as.character(.x)]][[as.character(.y)]]),
         S_valid = map(fold, ~ res_quict$S_valid[[.x]]),
         psi_measured = map(fold, ~ res_quict$psi_measured[[.x]]),
         psi_estimated = map2(fold, penalty,
                              ~ res_quict$psi_estimated[[as.character(.x)]][[as.character(.y)]]),
  )

# compute metrics
tib_quic <- tib_quic |>
  mutate(Rsquared = map2_dbl(psi_measured, psi_estimated,
                             ~ {
                               lm(as.numeric(.y) ~ as.numeric(.x)) |>
                                 summary() |>
                                 (\(x) x[["adj.r.squared"]])()
                             }),
         residuals = map2(psi_measured, psi_estimated,
                          ~ .y - .x),
         sum_abs_residuals = map_dbl(residuals, ~ sum(abs(.x))),
         FEV = map2(residuals, psi_measured, ~ frac_explained_var(.x, .y)),
         mean_FEV = map_dbl(FEV, ~ mean(.x)),
         loss_frobenius = map2_dbl(S_valid, S_train, ~loss_frob(.x, .y)),
         loss_quadratic = map2_dbl(S_valid, OM, ~loss_quad(.x, .y)))

# combine results of two previous runs with different rho
# tib_quic2 <- bind_rows(tib_quic,
#                        tib_quic1)
# qs::qsave(tib_quic2,
#           "data/intermediates/230920_cv/230921_tib_quic2.qs")

tib_quic <- qs::qread("data/intermediates/230920_cv/230921_tib_quic2.qs")

#~ Compile metrics ----

# Plot
tib_quic |>
  select(penalty, fold, Rsquared, sum_abs_residuals, mean_FEV, loss_frobenius, loss_quadratic) |>
  summarize(across(-fold, list(mean = mean, sd = sd)), .by = penalty) |>
  pivot_longer(-penalty,
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(metric = fct_inorder(metric)) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_point() + geom_line() + geom_errorbar() +
  scale_x_log10()

# export metrics as table
tib_quic |>
  select(penalty, fold, Rsquared, sum_abs_residuals, mean_FEV, loss_frobenius, loss_quadratic) |>
  summarize(across(-fold, mean), .by = penalty) |>
  arrange(desc(penalty)) |>
  mutate(Rsquared = round(Rsquared, 3),
         sum_abs_residuals = format(round(sum_abs_residuals), big.mark = ","),
         mean_FEV = round(mean_FEV, 3),
         loss_frobenius = round(loss_frobenius),
         loss_quadratic = round(loss_quadratic)) |> 
  clipr::write_clip(na = "NaN")







# huge (glasso) ----

rho_vals <- c(1, 0.81, 0.65, 0.53, 0.43, 0.34, 0.28, 0.22, 0.18, 0.15, 0.12, 0.1, 0.08, 0.06, 0.05) |>
  set_names()
fold_names <- sort(unique(folds)) |> set_names()



res_huge <- map(fold_names,
                \(fold){
                  
                  train_only <- mat_train[folds != fold, ] |>
                    huge::huge.npn(verbose = FALSE)
                  train_psi <- t(train_only[,1:nb_psi])
                  train_sf <- t(train_only[,(nb_psi+1):(nb_psi+nb_sf)])
                  
                  S_train <- cov(t(rbind(train_psi,train_sf)))
                  
                  # run estimation!
                  # # first pass: letting it choose its own lambda
                  # res_huge <- huge::huge(S_train,
                  #                        nlambda = 15,
                  #                        lambda.min.ratio = .05,
                  #                        method = "glasso",
                  #                        cov.output = TRUE)
                  
                  # second pass: forcing the lambda
                  res_huge <- huge::huge(S_train,
                                         lambda = rho_vals,
                                         method = "glasso",
                                         cov.output = TRUE)
                  
                  
                  # Extract results
                  rho_vals <- res_huge$lambda
                  f_OM <- res_huge$icov |>
                    lapply(\(OM) {dimnames(OM) <- dimnames(S_train); OM})
                  f_S_train <- res_huge$cov |>
                    lapply(\(S) {dimnames(S) <- dimnames(S_train); S})
                  
                  
                  
                  # get validation set
                  vald_only <- mat_train[folds == fold, ] |>
                    huge::huge.npn(verbose = FALSE)
                  vald_psi <- t(vald_only[,1:nb_psi])
                  vald_sf <- t(vald_only[,(nb_psi+1):(nb_psi+nb_sf)])
                  
                  S_valid <- cov(t(rbind(vald_psi,vald_sf)))
                  
                  
                  
                  f_psi_estimated <- lapply(f_OM,
                                            \(OM){
                                              OM21 <- OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
                                              OM11 <- OM[1:nb_psi, 1:nb_psi]
                                              
                                              W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                                              t(W) %*% vald_sf
                                            }) |>
                    set_names(rho_vals)
                  
                  
                  list(penalty = rho_vals,
                       S_train = f_S_train,
                       OM = f_OM,
                       S_valid = S_valid,
                       psi_measured = vald_psi,
                       psi_estimated = f_psi_estimated)
                },
                .progress = TRUE)

# qs::qsave(res_huge, "data/intermediates/230920_cv/230921_huge_fixed_lambda.qs")

# res_huge <- qs::qread("data/intermediates/230920_cv/230921_huge_fixed_lambda.qs")



#~ PSI prediction ----

par(mfrow = c(2,4), mar = c(2,2,2,2))
for(i in seq(1,15,2)) plot(res_huge$`1`$psi_measured,
                           res_huge$`1`$psi_estimated[[i]],
                           xlab = "",ylab = "",
                           main = paste("lambda =", res_huge$`1`$penalty[[i]] |> round(2)))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


# extract values
res_huget <- transpose(res_huge)



tib_huge <- res_huget |>
  as_tibble() |>
  unnest(c(penalty, S_train, OM, psi_estimated)) 

# compute metrics
tib_huge <- tib_huge |>
  mutate(Rsquared = map2_dbl(psi_measured, psi_estimated,
                             ~ {
                               lm(as.numeric(.y) ~ as.numeric(.x)) |>
                                 summary() |>
                                 (\(x) x[["adj.r.squared"]])()
                             }),
         residuals = map2(psi_measured, psi_estimated,
                          ~ .y - .x),
         sum_abs_residuals = map_dbl(residuals, ~ sum(abs(.x))),
         FEV = map2(residuals, psi_measured, ~ frac_explained_var(.x, .y)),
         mean_FEV = map_dbl(FEV, ~ mean(.x)),
         loss_frobenius = map2_dbl(S_valid, S_train, ~loss_frob(.x, .y)),
         loss_quadratic = map2_dbl(S_valid, OM, ~loss_quad(.x, .y)))


# qs::qsave(tib_huge,
#           "data/intermediates/230920_cv/230921_tib_huge.qs")

# tib_huge <- qs::qread("data/intermediates/230920_cv/230921_tib_huge.qs")



#~ Select penalties ----
# having run huge letting it select penalties, which are very similar from fold to fold,
# we can average them and use a consistent set of lambdas
plot(tib_huge$penalty, ylab = "penalty")
tibble(lambda = tib_huge$penalty,
       fold = rep(fold_names, each = 15),
       order = rep(1:15, times = length(fold_names))) |>
  summarize(lambda = mean(lambda),
            .by = order) |>
  pull(lambda) |>
  round(2) #|> clipr::write_clip()
#> not necessary to redo!



#~ Compile metrics ----
# Plot
tib_huge |>
  select(penalty, Rsquared, sum_abs_residuals, mean_FEV, loss_frobenius, loss_quadratic) |>
  summarize(across(everything(), list(mean = mean, sd = sd)), .by = penalty) |>
  pivot_longer(-penalty,
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(metric = fct_inorder(metric)) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_point() + geom_line() + geom_errorbar() +
  scale_x_log10()

# export metrics as table
tib_huge |>
  select(penalty, Rsquared, sum_abs_residuals, mean_FEV, loss_frobenius, loss_quadratic) |>
  summarize(across(everything(), mean), .by = penalty) |>
  arrange(desc(penalty)) |>
  mutate(Rsquared = round(Rsquared, 3),
         sum_abs_residuals = format(round(sum_abs_residuals), big.mark = ","),
         mean_FEV = round(mean_FEV, 3),
         loss_frobenius = round(loss_frobenius),
         loss_quadratic = round(loss_quadratic)) |> 
  clipr::write_clip(na = "NaN")










# huge (tiger) ----

rho_vals <- c(1, 0.85, 0.72, 0.61, 0.52, 0.44, 0.37, 0.32, 0.27, 0.23, 0.19, 0.16, 0.14, 0.12, 0.1) |>
  set_names()

fold_names <- sort(unique(folds)) |> set_names()



res_huge <- map(fold_names,
                \(fold){
                  
                  train_only <- mat_train[folds != fold, ] |>
                    huge::huge.npn(verbose = FALSE)
                  train_psi <- t(train_only[,1:nb_psi])
                  train_sf <- t(train_only[,(nb_psi+1):(nb_psi+nb_sf)])
                  
                  S_train <- cov(t(rbind(train_psi,train_sf)))
                  
                  # filter out rows/cols of 0
                  S_train2 <- S_train[rowSums(S_train) != 0, colSums(S_train) != 0]
                  
                  
                  # run estimation!
                  # # first pass: letting it choose its own lambda
                  # res_huge <- huge::huge(S_train2,
                  #                        nlambda = 15,
                  #                        method = "tiger")
                  # # second pass: forcing precomputed lambda
                  res_huge <- huge::huge(S_train2,
                                         lambda = rho_vals,
                                         method = "tiger")
                  
                  # Extract results
                  rho_vals <- res_huge$lambda
                  f_OM <- res_huge$icov |>
                    lapply(\(OM) {dimnames(OM) <- dimnames(S_train2); OM})
                  
                  
                  
                  # get validation set
                  vald_only <- mat_train[folds == fold, colnames(S_train2)] |>
                    huge::huge.npn(verbose = FALSE)
                  vald_psi <- t(vald_only[,startsWith(colnames(S_train2), "SE")])
                  vald_sf <- t(vald_only[,! startsWith(colnames(S_train2), "SE")])
                  
                  S_valid <- cov(t(rbind(vald_psi,vald_sf)))
                  
                  
                  
                  f_psi_estimated <- lapply(f_OM,
                                            \(OM){
                                              OM21 <- OM[! startsWith(colnames(S_train2), "SE"), startsWith(colnames(S_train2), "SE")]
                                              OM11 <- OM[startsWith(colnames(S_train2), "SE"), startsWith(colnames(S_train2), "SE")]
                                              
                                              W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                                              t(W) %*% vald_sf
                                            }) |>
                    set_names(rho_vals)
                  
                  
                  list(penalty = rho_vals,
                       OM = f_OM,
                       S_valid = S_valid,
                       psi_measured = vald_psi,
                       psi_estimated = f_psi_estimated)
                },
                .progress = TRUE)

# qs::qsave(res_huge, "data/intermediates/230920_cv/230922_huge_tiger_lambda.qs")

# res_huge <- qs::qread("data/intermediates/230920_cv/230922_huge_tiger_lambda.qs")



#~ PSI prediction ----

par(mfrow = c(2,4), mar = c(2,2,2,2))
for(i in seq(1,15,2)) plot(res_huge$`1`$psi_measured,
                           res_huge$`1`$psi_estimated[[i]],
                           xlab = "",ylab = "",
                           main = paste("lambda =", res_huge$`1`$penalty[[i]] |> round(2)))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


# extract values
res_huget <- transpose(res_huge)



tib_huge <- res_huget |>
  as_tibble() |>
  unnest(c(penalty, OM, psi_estimated)) 

# compute metrics
tib_huge <- tib_huge |>
  mutate(Rsquared = map2_dbl(psi_measured, psi_estimated,
                             ~ {
                               lm(as.numeric(.y) ~ as.numeric(.x)) |>
                                 summary() |>
                                 (\(x) x[["adj.r.squared"]])()
                             }),
         residuals = map2(psi_measured, psi_estimated,
                          ~ .y - .x),
         sum_abs_residuals = map_dbl(residuals, ~ sum(abs(.x))),
         FEV = map2(residuals, psi_measured, ~ frac_explained_var(.x, .y)),
         mean_FEV = map_dbl(FEV, ~ mean(.x)),
         loss_quadratic = map2_dbl(S_valid, OM, ~loss_quad(.x, .y)))


# qs::qsave(tib_huge,
#           "data/intermediates/230920_cv/230922_tib_huge_tiger.qs")

# tib_huge <- qs::qread("data/intermediates/230920_cv/230922_tib_huge_tiger.qs")



#~ Select penalties ----
# having run huge letting it select penalties, which are very similar from fold to fold,
# we can average them and use a consistent set of lambdas
plot(tib_huge$penalty, ylab = "penalty")
tibble(lambda = tib_huge$penalty,
       fold = rep(fold_names, each = 15),
       order = rep(1:15, times = length(fold_names))) |>
  summarize(lambda = mean(lambda),
            .by = order) |>
  pull(lambda) |>
  round(2) #|> clipr::write_clip()
#> not necessary to redo!



#~ Compile metrics ----
# Plot
tib_huge |>
  select(penalty, Rsquared, sum_abs_residuals, mean_FEV, loss_quadratic) |>
  summarize(across(everything(), list(mean = mean, sd = sd)), .by = penalty) |>
  pivot_longer(-penalty,
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(metric = fct_inorder(metric)) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_point() + geom_line() + geom_errorbar() +
  scale_x_log10()

# export metrics as table
tib_huge |>
  select(penalty, Rsquared, sum_abs_residuals, mean_FEV, loss_quadratic) |>
  summarize(across(everything(), mean), .by = penalty) |>
  arrange(desc(penalty)) |>
  mutate(Rsquared = round(Rsquared, 3),
         sum_abs_residuals = format(round(sum_abs_residuals), big.mark = ","),
         mean_FEV = round(mean_FEV, 3),
         loss_quadratic = round(loss_quadratic)) |> 
  clipr::write_clip(na = "NaN")







# ARACNE (make positive definite) ----


fold_names <- sort(unique(folds)) |> set_names()



res_arac <- map(fold_names,
                \(fold){
                  
                  train_only <- mat_train[folds != fold, ] |>
                    huge::huge.npn(verbose = FALSE)
                  train_psi <- t(train_only[,1:nb_psi])
                  train_sf <- t(train_only[,(nb_psi+1):(nb_psi+nb_sf)])
                  
                  S_train <- cov(t(rbind(train_psi,train_sf)))
                  
                  # run estimation!
                  mim <- minet::build.mim(t(rbind(train_psi,train_sf)))
                  mim_filt <- mim[rowSums(is.na(mim)) <= 1, colSums(is.na(mim)) <= 1]
                  arac <- minet::aracne(mim_filt)
                  arac2 <- corpcor::make.positive.definite(arac)
                  
                  
                  # get validation set
                  vald_only <- mat_train[folds == fold, colnames(arac2)] |>
                    huge::huge.npn(verbose = FALSE)
                  vald_psi <- t(vald_only[,startsWith(colnames(arac2), "SE")])
                  vald_sf <- t(vald_only[,! startsWith(colnames(arac2), "SE")])
                  
                  S_valid <- cov(t(rbind(vald_psi,vald_sf)))
                  
                  
                  
                  OM21 <- arac2[! startsWith(colnames(arac2), "SE"), startsWith(colnames(arac2), "SE")]
                  OM11 <- arac2[startsWith(colnames(arac2), "SE"), startsWith(colnames(arac2), "SE")]
                  
                  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                  pred <- t(W) %*% vald_sf
                  
                  
                  
                  list(S_train = S_train,
                       OM = arac2,
                       S_valid = S_valid,
                       psi_measured = vald_psi,
                       psi_estimated = pred)
                },
                .progress = TRUE)



#~ PSI prediction ----

plot(res_arac$`1`$psi_measured,
     res_arac$`1`$psi_estimated,
     xlab = "",ylab = "")


# extract values
res_aract <- transpose(res_arac)



tib_arac <- res_aract |>
  as_tibble()

# compute metrics
tib_arac <- tib_arac |>
  mutate(Rsquared = map2_dbl(psi_measured, psi_estimated,
                             ~ {
                               lm(as.numeric(.y) ~ as.numeric(.x)) |>
                                 summary() |>
                                 (\(x) x[["adj.r.squared"]])()
                             }),
         residuals = map2(psi_measured, psi_estimated,
                          ~ .y - .x),
         sum_abs_residuals = map_dbl(residuals, ~ sum(abs(.x))),
         FEV = map2(residuals, psi_measured, ~ frac_explained_var(.x, .y)),
         mean_FEV = map_dbl(FEV, ~ mean(.x)),
         loss_quadratic = map2_dbl(S_valid, OM, ~loss_quad(.x, .y)))


# qs::qsave(tib_arac,
#           "data/intermediates/230920_cv/230922_tib_arac.qs")







# export metrics as table
tib_arac |>
  select(Rsquared, sum_abs_residuals, mean_FEV, loss_quadratic) |>
  summarize(across(everything(), mean)) |>
  mutate(Rsquared = round(Rsquared, 3),
         sum_abs_residuals = format(round(sum_abs_residuals), big.mark = ","),
         mean_FEV = round(mean_FEV, 3),
         loss_quadratic = round(loss_quadratic)) |> 
  clipr::write_clip(na = "NaN")







# ARACNE (direct output) ----
# do not make positive definite (to keep matrix sparse),
# do not attempt all the metrics (will be used for ground truth only)

fold_names <- sort(unique(folds)) |> set_names()



res_arac <- map(fold_names,
                \(fold){
                  
                  train_only <- mat_train[folds != fold, ] |>
                    huge::huge.npn(verbose = FALSE)
                  train_psi <- t(train_only[,1:nb_psi])
                  train_sf <- t(train_only[,(nb_psi+1):(nb_psi+nb_sf)])
                  
                  S_train <- cov(t(rbind(train_psi,train_sf)))
                  
                  # run estimation!
                  mim <- minet::build.mim(t(rbind(train_psi,train_sf)))
                  mim_filt <- mim[rowSums(is.na(mim)) <= 1, colSums(is.na(mim)) <= 1]
                  arac <- minet::aracne(mim_filt)

                  
                  
                  
                  
                  list(S_train = S_train,
                       OM = arac)
                },
                .progress = TRUE)






# extract values
tib_arac <- transpose(res_arac) |>
  as_tibble()


# qs::qsave(tib_arac,
#           "data/intermediates/230920_cv/230922_tib_arac_raw.qs")











# DPM ----

# trying both with and without regularization (instead of penalties)

fold_names <- sort(unique(folds)) |> set_names()



res_dpm <- map(fold_names,
               \(fold){
                 
                 train_only <- mat_train[folds != fold, ] |>
                   huge::huge.npn(verbose = FALSE)
                 train_psi <- t(train_only[,1:nb_psi])
                 train_sf <- t(train_only[,(nb_psi+1):(nb_psi+nb_sf)])
                 
                 
                 ### Non-regularized ###
                 
                 # run estimation!
                 dpm_nonreg <- DPM::dpm(t(rbind(train_psi,train_sf)))
                 dimnames(dpm_nonreg) <- list(colnames(train_only),
                                              colnames(train_only))
                 dpm_nonreg <- dpm_nonreg[rowSums(is.na(dpm_nonreg)) <= 1, colSums(is.na(dpm_nonreg)) <= 1]
                 
                 
                 # get validation set
                 vald_only <- mat_train[folds == fold, colnames(dpm_nonreg)] |>
                   huge::huge.npn(verbose = FALSE)
                 vald_psi <- t(vald_only[, startsWith(colnames(dpm_nonreg), "SE")])
                 vald_sf <- t(vald_only[, !startsWith(colnames(dpm_nonreg), "SE")])
                 
                 S_valid_nonreg <- cov(t(rbind(vald_psi,vald_sf)))
                 
                 
                 OM21 <- dpm_nonreg[!startsWith(colnames(dpm_nonreg), "SE"),
                                    startsWith(colnames(dpm_nonreg), "SE")]
                 OM11 <- dpm_nonreg[startsWith(colnames(dpm_nonreg), "SE"),
                                    startsWith(colnames(dpm_nonreg), "SE")]
                 
                 W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                 pred_nonreg <- t(W) %*% vald_sf
                 
                 
                 
                 ### Regularized ###
                 # run estimation!
                 dpm_reg <- DPM::reg.dpm(t(rbind(train_psi,train_sf)))
                 dimnames(dpm_reg) <- list(colnames(train_only),
                                           colnames(train_only))
                 dpm_reg <- dpm_reg[rowSums(is.na(dpm_reg)) <= 1,
                                    colSums(is.na(dpm_reg)) <= 1]
                 
                 
                 # get validation set
                 vald_only <- mat_train[folds == fold, colnames(dpm_reg)] |>
                   huge::huge.npn(verbose = FALSE)
                 vald_psi <- t(vald_only[, startsWith(colnames(dpm_reg), "SE")])
                 vald_sf <- t(vald_only[, !startsWith(colnames(dpm_reg), "SE")])
                 
                 S_valid_reg <- cov(t(rbind(vald_psi,vald_sf)))
                 
                 
                 OM21 <- dpm_reg[!startsWith(colnames(dpm_reg), "SE"),
                                 startsWith(colnames(dpm_reg), "SE")]
                 OM11 <- dpm_reg[startsWith(colnames(dpm_reg), "SE"),
                                 startsWith(colnames(dpm_reg), "SE")]
                 
                 W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                 pred_reg <- t(W) %*% vald_sf
                 
                 
                 
                 list(OM = list(nonreg = dpm_nonreg, reg = dpm_reg),
                      S_valid = list(nonreg = S_valid_nonreg, reg = S_valid_reg),
                      psi_measured = vald_psi,
                      psi_estimated = list(nonreg = pred_nonreg, reg = pred_reg))
               },
               .progress = TRUE)

# qs::qsave(res_dpm, "data/intermediates/230920_cv/230922_dpm.qs")

# res_dpm <- qs::qread("data/intermediates/230920_cv/230922_dpm.qs")



#~ PSI prediction ----

par(mfrow = c(1,2))
plot(res_dpm$`1`$psi_measured,
     res_dpm$`1`$psi_estimated$nonreg,
     xlab = "",ylab = "")
plot(res_dpm$`1`$psi_measured,
     res_dpm$`1`$psi_estimated$reg,
     xlab = "",ylab = "")
par(mfrow = c(1,1))

# extract values
res_dpmt <- transpose(res_dpm)



tib_dpm <- res_dpmt |>
  as_tibble() |>
  mutate(fold = names(OM),
         .before = 1) |>
  unnest(c(OM, S_valid, psi_estimated)) |>
  mutate(penalty = names(OM),
         .before = 1)

# compute metrics
tib_dpm <- tib_dpm |>
  mutate(Rsquared = map2_dbl(psi_measured, psi_estimated,
                             ~{lm(as.numeric(.y) ~ as.numeric(.x)) |>
                                 summary() |>
                                 (\(x) x[["adj.r.squared"]])()
                             }),
         residuals = map2(psi_measured, psi_estimated,
                          ~ .y - .x),
         sum_abs_residuals = map_dbl(residuals, ~ sum(abs(.x))),
         FEV = map2(residuals, psi_measured, ~ frac_explained_var(.x, .y)),
         mean_FEV = map_dbl(FEV, ~ mean(.x)),
         loss_quadratic = map2_dbl(S_valid, OM, ~loss_quad(.x, .y)))


# qs::qsave(tib_dpm,
#           "data/intermediates/230920_cv/230922_tib_dpm.qs")

# tib_huge <- qs::qread("data/intermediates/230920_cv/230921_tib_huge.qs")




#~ Compile metrics ----

# Plot
tib_dpm |>
  select(penalty, Rsquared, sum_abs_residuals, mean_FEV, loss_quadratic) |>
  summarize(across(everything(), list(mean = mean, sd = sd)), .by = penalty) |>
  pivot_longer(-penalty,
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(metric = fct_inorder(metric)) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_point() + geom_errorbar()


# export metrics as table
tib_dpm |>
  select(penalty, Rsquared, sum_abs_residuals, mean_FEV, loss_quadratic) |>
  summarize(across(everything(), mean), .by = penalty) |>
  arrange(desc(penalty)) |>
  mutate(Rsquared = round(Rsquared, 3),
         sum_abs_residuals = format(round(sum_abs_residuals), big.mark = ","),
         mean_FEV = round(mean_FEV, 3),
         loss_quadratic = round(loss_quadratic)) |> 
  clipr::write_clip(na = "NaN")












# scio ----

rho_vals <- c(10, 5, 2, 1, .5, .1) |> set_names()
rho_vals <- c(.05, .025, .01, .0075) |> set_names()


fold_names <- sort(unique(folds)) |> set_names()



res_scio <- map(fold_names,
                \(fold){
                  
                  train_only <- mat_train[folds != fold, ] |>
                    huge::huge.npn(verbose = FALSE)
                  train_psi <- t(train_only[,1:nb_psi])
                  train_sf <- t(train_only[,(nb_psi+1):(nb_psi+nb_sf)])
                  
                  S_train <- cov(t(rbind(train_psi,train_sf)))
                  
                  # filter out rows/cols with no variance (that make OM non-invertible)
                  S_train <- S_train[matrixStats::rowVars(S_train) > 0,
                                     matrixStats::colVars(S_train) > 0]
                  
                  is_psi <- startsWith(colnames(S_train), "SE_")
                  is_sf <- !startsWith(colnames(S_train), "SE_")
                  
                  
                  # run estimation!
                  f_OM <- map(rho_vals,
                              ~ {
                                res <- scio::scio(S_train, lambda = .x)
                                res <- res[["w"]]
                                dimnames(res) <- dimnames(S_train)
                                res}
                  )
                  
                  
                  
                  # Extract results
                  f_S_train <- map(f_OM, possibly(\(x) solve(x),
                                                  otherwise = NA))
                  
                  # get validation set
                  vald_only <- mat_train[folds == fold, colnames(S_train)] |>
                    huge::huge.npn(verbose = FALSE)
                  vald_psi <- t(vald_only[, is_psi])
                  vald_sf <- t(vald_only[, is_sf])
                  
                  S_valid <- cov(t(rbind(vald_psi,vald_sf)))
                  
                  
                  f_psi_estimated <- map(f_OM,
                                         \(OM){
                                           OM21 <- OM[is_sf, is_psi]
                                           OM11 <- OM[is_psi, is_psi]
                                           
                                           W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                                           t(W) %*% vald_sf
                                         })
                  
                  
                  list(S_train = f_S_train,
                       OM = f_OM,
                       S_valid = S_valid,
                       psi_measured = vald_psi,
                       psi_estimated = f_psi_estimated)
                },
                .progress = TRUE)

# qs::qsave(res_scio, "data/intermediates/230920_cv/230925_scio2.qs")







#~ PSI prediction ----

par(mfrow = c(2,3), mar = c(2,2,2,2))
for(i in seq_along(res_scio$`1`$OM)) plot(res_scio$`1`$psi_measured,
                                          res_scio$`1`$psi_estimated[[i]], 
                                          xlab = "",ylab = "", main = paste("rho =",rho_vals[[i]]))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)


# extract values
res_sciot <- transpose(res_scio)

tib_scio <- expand_grid(penalty = rho_vals, fold = fold_names) |>
  mutate(S_train = map2(fold, penalty,
                        ~ res_sciot$S_train[[as.character(.x)]][[as.character(.y)]]),
         OM = map2(fold, penalty,
                   ~ res_sciot$OM[[as.character(.x)]][[as.character(.y)]]),
         S_valid = map(fold, ~ res_sciot$S_valid[[.x]]),
         psi_measured = map(fold, ~ res_sciot$psi_measured[[.x]]),
         psi_estimated = map2(fold, penalty,
                              ~ res_sciot$psi_estimated[[as.character(.x)]][[as.character(.y)]]),
  )

# compute metrics
tib_scio <- tib_scio |>
  mutate(Rsquared = map2_dbl(psi_measured, psi_estimated,
                             ~ {
                               lm(as.numeric(.y) ~ as.numeric(.x)) |>
                                 summary() |>
                                 (\(x) x[["adj.r.squared"]])()
                             }),
         residuals = map2(psi_measured, psi_estimated,
                          ~ .y - .x),
         sum_abs_residuals = map_dbl(residuals, ~ sum(abs(.x))),
         FEV = map2(residuals, psi_measured, ~ frac_explained_var(.x, .y)),
         mean_FEV = map_dbl(FEV, ~ mean(.x)),
         loss_frobenius = map2_dbl(S_valid, S_train, ~loss_frob(.x, .y)),
         loss_quadratic = map2_dbl(S_valid, OM, ~loss_quad(.x, .y)))


# qs::qsave(tib_scio,
#           "data/intermediates/230920_cv/230925_tib_scio2.qs")
# tib_scio <- qs::qread("data/intermediates/230920_cv/230925_tib_scio2.qs")
# tib_scio1 <- qs::qread("data/intermediates/230920_cv/230925_tib_scio.qs")

# tib3 <- bind_rows(tib_scio, tib_scio1)
# qs::qsave(tib3, "data/intermediates/230920_cv/230925_tib_scio_combined.qs")


# tib_scio <- qs::qread("data/intermediates/230920_cv/230925_tib_scio_combined.qs")


#~ Compile metrics ----

# Plot
tib_scio |>
  select(penalty, fold, Rsquared, sum_abs_residuals, mean_FEV, loss_frobenius, loss_quadratic) |>
  summarize(across(-fold, list(mean = mean, sd = sd)), .by = penalty) |>
  pivot_longer(-penalty,
               names_to = c("metric", "type"),
               names_pattern = "(.+)_(mean|sd)$",
               values_to = "value") |>
  pivot_wider(names_from = "type",
              values_from = "value") |>
  mutate(metric = fct_inorder(metric)) |>
  ggplot(aes(x = penalty, y = mean, ymin = mean - sd, ymax = mean + sd)) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y") +
  geom_point() + geom_line() + geom_errorbar() +
  scale_x_log10()

# export metrics as table
tib_scio |>
  select(penalty, fold, Rsquared, sum_abs_residuals, mean_FEV, loss_frobenius, loss_quadratic) |>
  summarize(across(-fold, mean), .by = penalty) |>
  arrange(desc(penalty)) |>
  mutate(Rsquared = round(Rsquared, 3),
         sum_abs_residuals = format(round(sum_abs_residuals), big.mark = ","),
         mean_FEV = round(mean_FEV, 3),
         loss_frobenius = round(loss_frobenius),
         loss_quadratic = round(loss_quadratic)) |> 
  clipr::write_clip(na = "NaN")






























