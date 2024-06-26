# copied from graph_power4, run on cluster to use many cores


# Based on graph_power_cv and graph_power_cv3; only NPN/imputation within each fold, test methods


message("Starting, ", date())

# Inits ----

suppressPackageStartupMessages(library(tidyverse))
library(furrr)

plan(multicore, workers = 6)


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
  select(event_id, sample_id, PSI) |>
  pivot_wider(id_cols = sample_id,
              names_from = event_id,
              values_from = PSI,
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
    separate_wider_delim(event_id_percount,
                         delim = ".",
                         names = c("event_id", NA)) |>
    group_by(sf_id, event_id) |>
    summarize(coefficient = max(coefficient),
              .groups = "drop") |>
    left_join(all_interactions_by_event,
              by = c("event_id", "sf_id" = "sf_tx")) |>
    mutate(literature = replace_na(literature, FALSE))
}




# ********** ----










# QUIC ----


rho_vals <- c(10, 5, 2, 1, .5, .1) |>
  set_names()

fold_names <- sort(unique(folds)) |> set_names()

message("---- Starting!!")

#~ prepare data -----
res_quic1 <- expand_grid(fold = fold_names,
                        permutation = 0:50) |>
  mutate(
    # training set
    psi_train = future_map2(fold, permutation,
                     ~{
                       out <- mat_train[folds != .x, 1:nb_psi] |>
                         huge::huge.npn(verbose = FALSE)
                       if(.y){
                         rownm <- rownames(out)
                         out <- apply(out, 2, sample)
                         rownames(out) <- rownm
                       }
                       out
                     },
		    .options = furrr_options(seed = TRUE)
		    ),
    sf_train = future_map(fold,
                   ~{
                     mat_train[folds != .x, (nb_psi+1):(nb_psi+nb_sf)] |>
                       huge::huge.npn(verbose = FALSE)
                   }),
    S_train = future_map2(psi_train, sf_train,
                   ~ cov(cbind(.x,.y))),
    #validation set
    psi_valid = future_map(fold,
                    ~{
                      mat_train[folds == .x, 1:nb_psi] |>
                        huge::huge.npn(verbose = FALSE)
                    }),
    sf_valid = future_map(fold,
                   ~{
                     mat_train[folds == .x, (nb_psi+1):(nb_psi+nb_sf)] |>
                       huge::huge.npn(verbose = FALSE)
                   }),
    
    S_valid = future_map2(psi_valid, sf_valid,
                   ~ cov(cbind(.x,.y)))
  ) 
  #~ estimate precision matrix! -----
res_quic2 <- res_quic1 |> mutate(fit = future_map(S_train,
                        ~ QUIC::QUIC(.x,
                                     rho = 1, path = rho_vals,
                                     msg = 0),
                        .progress = TRUE,
			.options = furrr_options(seed = TRUE)))
                        
  # extract estimates
res_quic3 <- res_quic2  |> mutate(OM = future_map2(fit, S_train,
                   \(.fit, .S_train){
                     map(seq_along(rho_vals) |> set_names(rho_vals),
                         ~ {
                           OM <- .fit[["X"]][,, .x ]
                           dimnames(OM) <- dimnames(.S_train)
                           OM
                         })
                   }),
         S_train_hat = future_map2(fit, S_train,
                            \(.fit, .S_train){
                              map(seq_along(rho_vals) |> set_names(rho_vals),
                                  ~ {
                                    OM <- .fit[["W"]][,, .x ]
                                    dimnames(OM) <- dimnames(.S_train)
                                    OM
                                  })
                            }))
res_quic <- res_quic3 |>
  unnest(c(OM, S_train_hat))|>
  mutate(penalty = names(OM) |> as.numeric(),
         .before = 2)

#~ compute CV estimate of psi ----
tib_quic <- res_quic |>
  mutate(psi_estimated = future_map2(OM, sf_valid,
                              \(.OM, .sf_valid){
                                OM21 <- .OM[(nb_psi + 1):(nb_psi + nb_sf), 1:nb_psi]
                                OM11 <- .OM[1:nb_psi, 1:nb_psi]
                                
                                W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
                                t(t(W) %*% t(.sf_valid))
                              })) |>
  # compute metrics
  mutate(Rsquared = future_map2_dbl(psi_valid, psi_estimated,
                             ~ {
                               lm(as.numeric(.y) ~ as.numeric(.x)) |>
                                 summary() |>
                                 (\(x) x[["adj.r.squared"]])()
                             }),
         residuals = future_map2(psi_valid, psi_estimated,
                          ~ .y - .x),
         sum_abs_residuals = future_map_dbl(residuals, ~ sum(abs(.x))),
         FEV = future_map2(residuals, psi_valid, ~ frac_explained_var(.x, .y)),
         mean_FEV = future_map_dbl(FEV, ~ mean(.x)),
         loss_frobenius = future_map2_dbl(S_valid, S_train_hat, ~loss_frob(.x, .y)),
         loss_quadratic = future_map2_dbl(S_valid, OM, ~loss_quad(.x, .y)),
         # process ground truth
         adj = future_map(OM, get_coefs_from_OM,
                   .progress = TRUE),
         prop_non_zero_coefs_litt = future_map_dbl(adj,
                                            ~ mean(.x$coefficient[.x$literature] != 0)),
         prop_non_zero_coefs_nonlitt = future_map_dbl(adj,
                                               ~ mean(.x$coefficient[! .x$literature] != 0)))

message("Saving, ", date())

qs::qsave(tib_quic, file.path(outdir,
                              "231107_quic_50perm.qs"))


message("All done, ", date())

