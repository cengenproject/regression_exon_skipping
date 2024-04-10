# copied from graph_power4, run on cluster to use many cores


# Based on graph_power_cv and graph_power_cv3; only NPN/imputation within each fold, test methods


message("Starting, ", date())

# Inits ----

suppressPackageStartupMessages(library(tidyverse))
library(getopt)



if(! interactive()){
  spec <- matrix(c(
    'date', 'd', 1, 'character',
    'exonsInput', 'e', 1, 'character',
    'transformation', 't', 1, 'character',
    'imputation', 'i', 1, 'character',
    'permutations', 'p', 1, 'character',
    'penalties', 'n', 1, 'character',
    'algo', 'a', 1, 'character',
    'knn_k', 'k', 2, 'integer',
    'save', 's', 0, 'logical'
  ), byrow=TRUE, ncol=4)
  
  params <- getopt(spec)
  
} else{
  # Options for interactive
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
}

# check and parse arguments
if(str_detect(params$permutations, "^[0-9]$")){
  params$permutations <- as.integer(params$permutations)
} else if(str_detect(params$permutations, "^[0-9]+:[0-9]+$")){
  perms <- str_split_1(params$permutations, ":")
  params$permutations <- seq(from = perms[[1]], to = perms[[2]])
  stopifnot(is.integer(params$permutations) && length(params$permutations) > 1)
} else{
  stop("Failed to parse 'permutations' parameter")
}


stopifnot(str_detect(params$penalties, "^[0-9c()., ]+$"))
params$penalties <- eval(str2expression(params$penalties))
stopifnot(is.numeric(params$penalties) && length(params$penalties) > 0)


params$exonsInput <- match.arg(params$exonsInput,
                               choices = c("PSI", "counts"))

params$transformation <- match.arg(params$transformation,
                                   choices = c("npnshrink", "npntrunc", "zscore"))

params$imputation <- match.arg(params$imputation,
                               choices = c("median", "knn"))

params$algo <- match.arg(params$algo,
                         choices = c("QUIC", "glasso", "CLIME", "SCIO"))

if(params$imputation != "knn"){
  if(! is.null(params$knn_k)){
    warning("knn_k will be ignored as imputation method is ", params$imputation)
  }
} else{
  if(is.null(params$knn_k)) params$knn_k <- 10
}




#~ Load ----

source("R/loss_functions.R")
source("R/functions_steps.R")

outdir <- "data/graph_power4/outputs"


datadir <- "data/graph_power4/inputs/240410_precomputed/"

nb_se <- qs::qread(file.path(datadir, "nb_se.qs"))
nb_sf <- qs::qread(file.path(datadir, "nb_sf.qs"))
mat_interactions_lit <- qs::qread(file.path(datadir, "mat_interactions_lit.qs"))


# We're not doing CV, but reusing CV scripts: as if test was the second fold
mat_all_train <- switch(params$exonsInput,
                    PSI = qs::qread(file.path(datadir, "mat_train_psi.qs")),
                    counts = qs::qread(file.path(datadir, "mat_train_cnt.qs")))
mat_all_test <- switch(params$exonsInput,
                        PSI = qs::qread(file.path(datadir, "mat_test_psi.qs")),
                        counts = qs::qread(file.path(datadir, "mat_test_cnt.qs")))

mat_train <- rbind(mat_all_train, mat_all_test)
folds <- rep(c(1,2), times = c(nrow(mat_all_train), nrow(mat_all_test)))


nb_se <- switch (params$exonsInput,
                 PSI = nb_se,
                 counts = 2*nb_se
)






# Analysis ----

rho_vals <- params$penalties |>
  set_names()

# single fold: `1` is the training data, `2` the test data
fold_names <- 2 |> set_names()

message("---- Starting!!")

set.seed(sum(params$permutations))
#~ prepare data -----
res_quic1 <- expand_grid(fold = fold_names,
                         permutation = params$permutations)

res_quic1$se_train_t <- map2(res_quic1$fold, res_quic1$permutation,
                             extract_transform_se_train)

res_quic1$sf_train_t <- map(res_quic1$fold,
                            extract_transform_sf_train)

res_quic1$S_train_t <- map2(res_quic1$se_train_t, res_quic1$sf_train_t,
                            compute_S)

res_quic1$se_valid_u = map(res_quic1$fold, 
                           extract_se_valid)


res_quic1$psi_valid_u <- switch(
  params$exonsInput,
  PSI = res_quic1$se_valid_u,
  counts = map(res_quic1$se_valid_u, reconstruct_psi_from_counts)
)


res_quic1$se_valid_t = map2(res_quic1$se_valid_u, res_quic1$se_train_t,
                            transform_from_prev)

res_quic1$sf_valid_t <- map2(res_quic1$fold, res_quic1$sf_train_t,
                             extract_transform_sf_valid)

res_quic1$S_valid_t <- map2(res_quic1$se_valid_t, res_quic1$sf_valid_t,
                            compute_S)


#~ estimate precision matrix! -----
message("  Estimate precision matrix")
res_quic1$fit <- map(res_quic1$S_train_t,
                     estimate_precision_mat,
                     .progress = TRUE)


# extract estimates
message("  Extract estimates")
res_quic1$OM_train <- map2(res_quic1$fit, res_quic1$S_train_t,
                           extract_precision_mat_estimate)

res_quic1$S_train_hat_t <- map2(res_quic1$fit, res_quic1$S_train_t,
                                extract_S_train_estimate)



res_quic <- res_quic1 |>
  unnest(c(OM_train, S_train_hat_t))|>
  mutate(penalty = names(OM_train) |> as.numeric(),
         .before = 2)


#~ compute results ----
message("  Compute results")

#~~ adjacency ----

res_quic$adj <- switch(
  params$exonsInput,
  
  PSI = map(res_quic$OM_train, extract_adj_mat),
  
  counts = res_quic$OM_train |>
    map(extract_adj_mat) |>
    map(reconstruct_adj_psi)
)




#~~ PSI hat ----

res_quic$se_valid_hat_t <- map2(res_quic$OM_train,
                                res_quic$sf_valid_t,
                                estimate_se)

res_quic$se_valid_hat_u <- map2(res_quic$se_valid_hat_t,
                                res_quic$se_train_t,
                                untransform_se_hat)

res_quic$psi_valid_hat_u <- switch(
  params$exonsInput,
  
  PSI = res_quic$se_valid_hat_u,
  counts = map(res_quic$se_valid_hat_u,
               reconstruct_psi_from_counts)
)









#~ compute metrics ----
message("  Compute metrics")

#~~ fitting loss ----
res_quic$loss_frobenius = map2_dbl(res_quic$S_valid_t,
                                   res_quic$S_train_hat_t,
                                   loss_frob)
res_quic$loss_quadratic = map2_dbl(res_quic$S_valid_t,
                                   res_quic$OM_train,
                                   loss_quad)

res_quic$bias_loss_frobenius = map2_dbl(res_quic$S_train_t,
                                        res_quic$S_train_hat_t,
                                        loss_frob)
res_quic$bias_loss_quadratic = map2_dbl(res_quic$S_train_t,
                                        res_quic$OM_train,
                                        loss_quad)


#~~ reconstruction ----
res_quic$Rsquared <- map2_dbl(res_quic$psi_valid_u,
                              res_quic$psi_valid_hat_u,
                              compute_metric_rsquared)


res_quic$mean_FEV <- map2_dbl(res_quic$psi_valid_u,
                              res_quic$psi_valid_hat_u,
                              compute_metric_FEV)





#~~ bio relevance ----

res_quic$literature_TPR = map_dbl(res_quic$adj,
                                  compute_TPR)
res_quic$literature_FPR = map_dbl(res_quic$adj,
                                  compute_FPR)

res_quic$power_law <- map_dbl(res_quic$adj,
                              mat_power_law)

res_quic$sparsity <- map_dbl(res_quic$adj, adj_sparsity)



# Save ----

if(params$imputation == "knn"){
  k_str <- paste0(params$imputation, "_", "k", params$knn_k)
} else{
  k_str <- params$imputation
}

out_name <- paste(params$date, params$algo, params$exonsInput,
                  params$transformation, k_str,
                  length(params$permutations), length(params$penalties),
                  sep = "_")


message("Saving as, ", out_name, " at ", date())

if(! is.null(params$save)){
  qs::qsave(res_quic, file.path(outdir,
                                paste0(out_name, ".qs")))
}


res_quic |>
  dplyr::select(penalty, fold, permutation, Rsquared,
                mean_FEV, loss_frobenius, loss_quadratic,
                bias_loss_frobenius, bias_loss_quadratic,
                literature_TPR, literature_FPR, power_law, sparsity) |>
  readr::write_csv(file.path(outdir,
                             paste0(out_name, ".csv")))


message("All done, ", date())

