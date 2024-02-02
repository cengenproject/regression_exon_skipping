# copied from graph_power4, run on cluster to use many cores


# Based on graph_power_cv and graph_power_cv3; only NPN/imputation within each fold, test methods


message("Starting, ", date())

# Inits ----

suppressPackageStartupMessages(library(tidyverse))


library(wbData)


gids <- wb_load_gene_ids(281)
tx2g <- wb_load_tx2gene(281)

source("R/loss_functions.R")

datadir <- "data/graph_power4/inputs"
outdir <- "data/graph_power4/outputs"


#~ Load data ----
quantifs_filtered <- qs::qread(file.path(datadir, "230206_preprocessed_quantifs_filtered.qs"))
sf_expression <- qs::qread(file.path(datadir, "230206_preprocessed_sf_expression.qs")) |>
  filter(transcript_id != "R07E5.14.2")


#~ Prepare data ----
message("---- Prepare data")

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


# Get Ground truth ----
message("---- Get ground truth")


events_coords <- read_tsv(file.path(datadir, "221111_events_coordinates.tsv"),
                          show_col_types = FALSE) |>
  select(event_id, gene_id)

all_interactions <- read_tsv(file.path(datadir, "sf_targets_v3.tsv"),
                                    show_col_types = FALSE)

all_interactions_by_event <- all_interactions |>
  select(SF, targets) |>
  distinct() |>
  left_join(events_coords,
            by = join_by(targets == gene_id),
            relationship = "many-to-many") |>
  filter(!is.na(event_id)) |>
  select(event_id, SF) |>
  mutate(sf_tx = wb_g2tx(SF, tx2g)) |>
  select(-SF) |>
  unnest(sf_tx) |>
  distinct() |>
  add_column(literature = TRUE)




get_coefs_from_OM <- function(OM){
  
  OM[startsWith(rownames(OM), "SE_"),
     !startsWith(colnames(OM), "SE_")] |>
    as.data.frame() |>
    rownames_to_column("event_id_percount") |>
    pivot_longer(cols = -event_id_percount,
                 names_to = "sf_id",
                 values_to = "coefficient") |>
    mutate(event_id = str_extract(event_id_percount, "^SE_[0-9]+")) |>
    select(-event_id_percount) |>
    group_by(sf_id, event_id) |>
    summarize(coefficient = max(coefficient),
              .groups = "drop") |>
    left_join(all_interactions_by_event,
              by = c("event_id", "sf_id" = "sf_tx")) |>
    mutate(literature = replace_na(literature, FALSE))
}

recompute_psi_from_counts <- function(mat_N){
  stopifnot(identical(str_split_i(colnames(mat_N)[c(TRUE,FALSE)], "\\.", 1),
                      str_split_i(colnames(mat_N)[c(FALSE,TRUE)], "\\.", 1)))
  
  stopifnot(identical(str_replace(colnames(mat_N)[c(TRUE,FALSE)],
                                  "\\.Nincl", "\\.Nexcl"),
                      colnames(mat_N)[c(FALSE,TRUE)]))
  
  mat_Nincl <- mat_N[,c(TRUE,FALSE)]
  colnames(mat_Nincl) <- str_remove(colnames(mat_Nincl), "\\.Nincl$")
  mat_Nexcl <- mat_N[,c(FALSE,TRUE)]
  colnames(mat_Nexcl) <- str_remove(colnames(mat_Nexcl), "\\.Nexcl$")
  
  all.equal(dimnames(mat_Nincl), dimnames(mat_Nexcl))
  
  
  mat_Nincl/(mat_Nincl + mat_Nexcl)
}


# ********** ----










# QUIC ----


rho_vals <- c(10, 5, 2, 1, .5, .1, .05) |>
  set_names()


fold_names <- sort(unique(folds)) |> set_names()

message("---- Starting!!")

#~ prepare data -----
res_quic1 <- expand_grid(fold = fold_names,
                         permutation = 0)

res_quic1$psi_train_t <- map2(res_quic1$fold, res_quic1$permutation,
                              ~{
                                out <- mat_train[folds != .x, 1:nb_psi] |>
                                  projectNPN::transform_npn_shrinkage()
                                if(.y){
                                  rownm <- rownames(out$mat)
                                  out$mat <- apply(out$mat, 2, sample)
                                  rownames(out$mat) <- rownm
                                }
                                out
                              }
)

res_quic1$sf_train_t <- map(res_quic1$fold,
                            ~{
                              mat_train[folds != .x, (nb_psi+1):(nb_psi+nb_sf)] |>
                                projectNPN::transform_npn_shrinkage()
                            })

res_quic1$S_train_t <- map2(res_quic1$psi_train_t, res_quic1$sf_train_t,
                            ~ cov(cbind(.x$mat, .y$mat)))

res_quic1$psi_valid_u = map(res_quic1$fold, 
                            ~{
                              mat_train[folds == .x, 1:nb_psi]
                            })

res_quic1$psi_valid_t = map2(res_quic1$psi_valid_u, res_quic1$psi_train_t,
                             ~{
                               projectNPN::transform_npn_shrinkage(.x, .y[["parameters"]])
                             })

res_quic1$sf_valid_t <- map2(res_quic1$fold, res_quic1$sf_train_t,
                             ~{
                               mat_train[folds == .x, (nb_psi+1):(nb_psi+nb_sf)] |>
                                 projectNPN::transform_npn_shrinkage(.y[["parameters"]])
                             })

res_quic1$S_valid_t <- map2(res_quic1$psi_valid_t, res_quic1$sf_valid_t,
                            ~ cov(cbind(.x$mat, .y$mat)))


#~ estimate precision matrix! -----

res_quic1$fit <- map(res_quic1$S_train_t,
                        ~ QUIC::QUIC(.x,
                                     rho = 1, path = rho_vals,
                                     msg = 0),
                        .progress = TRUE)

message("Done estimating precision matrix")

# extract estimates
res_quic1$OM_train <- map2(res_quic1$fit, res_quic1$S_train_t,
                           \(.fit, .S_train){
                             map(seq_along(rho_vals) |> set_names(rho_vals),
                                 ~ {
                                   OM <- .fit[["X"]][,, .x ]
                                   dimnames(OM) <- dimnames(.S_train)
                                   OM
                                 })
                           })

res_quic1$S_train_hat_t <- map2(res_quic1$fit, res_quic1$S_train_t,
                                \(.fit, .S_train){
                                  map(seq_along(rho_vals) |> set_names(rho_vals),
                                      ~ {
                                        S <- .fit[["W"]][,, .x ]
                                        dimnames(S) <- dimnames(.S_train)
                                        S
                                      })
                                })

message("Done extracting fits")

res_quic <- res_quic1 |>
  unnest(c(OM_train, S_train_hat_t))|>
  mutate(penalty = names(OM_train) |> as.numeric(),
         .before = 2)

#~ compute CV estimate of psi ----
res_quic$psi_valid_hat_t <- map2(res_quic$OM_train, res_quic$sf_valid_t,
                                 \(.OM, .sf_valid){
                                   OM21 <- .OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
                                   OM11 <- .OM[1:nb_psi, 1:nb_psi]
                                   
                                   W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                                   t(t(W) %*% t(.sf_valid$mat))
                                 })

res_quic$psi_valid_hat_u <- map2(res_quic$psi_valid_hat_t,
                                 res_quic$psi_train_t,
                                 ~{
                                   projectNPN::rev_npn_shrinkage(.x,
                                                                 .y$parameters)
                                 })



#~ Convert back to PSI ----


res_quic$rpsi_valid_hat_u <- map(res_quic$psi_valid_hat_u,
                                 recompute_psi_from_counts)

res_quic$rpsi_valid_u <- map(res_quic$psi_valid_u,
                                 recompute_psi_from_counts)



#~ compute metrics ----
res_quic$Rsquared <- map2_dbl(res_quic$rpsi_valid_u, res_quic$rpsi_valid_hat_u,
                             ~ {
                               lm(as.numeric(.y) ~ as.numeric(.x)) |>
                                 summary() |>
                                 (\(x) x[["adj.r.squared"]])()
                             })
res_quic$residuals = map2(res_quic$rpsi_valid_u, res_quic$rpsi_valid_hat_u,
                          ~ .y - .x)
res_quic$sum_abs_residuals = map_dbl(res_quic$residuals, ~ sum(abs(.x), na.rm = TRUE))
res_quic$FEV = map2(res_quic$residuals, res_quic$rpsi_valid_u, ~ frac_explained_var(.x, .y, na.rm = TRUE))
res_quic$mean_FEV = map_dbl(res_quic$FEV, ~ mean(.x))

res_quic$loss_frobenius = map2_dbl(res_quic$S_valid_t, res_quic$S_train_hat_t, ~loss_frob(.x, .y))
res_quic$loss_quadratic = map2_dbl(res_quic$S_valid_t, res_quic$OM_train, ~loss_quad(.x, .y))

# process ground truth
res_quic$adj = map(res_quic$OM_train, get_coefs_from_OM,
                   .progress = TRUE)
res_quic$prop_non_zero_coefs_litt = map_dbl(res_quic$adj,
                                            ~ mean(.x$coefficient[.x$literature] != 0))
res_quic$prop_non_zero_coefs_nonlitt = map_dbl(res_quic$adj,
                                               ~ mean(.x$coefficient[! .x$literature] != 0))

# Save ----

out_name <- "240202_recompandrevertpsi_noperm_7penalties"

message("Saving as, ", out_name, " at ", date())

qs::qsave(res_quic, file.path(outdir,
                              paste0(out_name, ".qs")))

res_quic |>
  dplyr::select(penalty, fold, permutation, Rsquared, sum_abs_residuals,
                mean_FEV, loss_frobenius, loss_quadratic,
                prop_non_zero_coefs_litt, prop_non_zero_coefs_nonlitt) |>
  readr::write_csv(file.path(outdir,
                              paste0(out_name, ".csv")))


message("All done, ", date())

