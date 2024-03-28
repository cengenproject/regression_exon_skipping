# steps that get called from main script
# Note they make heavy use of global variables!


# Generic functions depending on params

impute <- switch(params$imputation,
                 median = impute_median,
                 knn = imute_knn)


transform_fwd <- switch(params$transformation,
                        npnshrink = projectNPN::transform_npn_shrinkage,
                        npntrunc = projectNPN::transform_npn_truncation,
                        zscore = projectNPN::transform_zscore)
  

transform_rev <- switch(params$transformation,
                        npnshrink = projectNPN::reverse_npn_shrinkage,
                        npntrunc = projectNPN::reverse_npn_truncation,
                        zscore = projectNPN::reverse_transform_zscore)





# Each step

extract_transform_se_train <- function(.fold, .permutation){
  
  out <- mat_train[folds != .fold, 1:nb_se] |>
    impute() |>
    transform_fwd()
  
  if(.permutation){
    rownm <- rownames(out$mat)
    out$mat <- apply(out$mat, 2, sample)
    rownames(out$mat) <- rownm
  }
  
  out
}



extract_transform_sf_train <- function(.fold){
  
  mat_train[folds != .fold, (nb_se+1):(nb_se+nb_sf)] |>
    transform_fwd(na = "center")
}


compute_S <- function(se, sf){
  
  cov(cbind(se$mat, sf$mat))
}




extract_se_valid <- function(.fold){
  
  mat_train[folds == .fold, 1:nb_se]
}



transform_from_prev <- function(untransformed, prev_transformed){
  
  untransformed |>
    transform_fwd(prev_transformed[["parameters"]], na = "keep")
}



extract_transform_sf_valid <- function(.fold, prev_transformed){
  
  mat_train[folds == .fold, (nb_se+1):(nb_se+nb_sf)] |>
    transform_fwd(prev_transformed[["parameters"]], na = "center")
}



estimate_precision_mat <- switch(
  params$algo,
  
  QUIC = function(.S){
    QUIC::QUIC(.S,
               rho = 1, path = rho_vals,
               msg = 0)
  },
  
  glasso = function(.S){
    map(rho_vals, ~glasso::glasso(.S, rho = .x))
  },
  
  CLIME = function(.S){
   flare::sugm(.S,
               lambda = rho_vals,
               method = "clime",
               verbose = FALSE)
  },
  
  SCIO = function(.S){
    map(rho_vals, ~scio::scio(.S, lambda = .x))
  }
)


# Notes on flare output:
# sigma is the input matrix, sigma2 has "perturb" added on the diagonal
# icov is the symmetrized version of icov1 (specifically, for any element 
# outside of diag, icov keeps the smallest of the sym, checking the source)


extract_precision_mat_estimate <- switch(
  params$algo,
  QUIC = function(.fit, .S_train){
    
    map(seq_along(rho_vals) |> set_names(rho_vals),
        ~ {
          OM <- .fit[["X"]][,, .x ]
          dimnames(OM) <- dimnames(.S_train)
          OM
        })
  },
  glasso = function(.fit, .S_train){
    
    map(seq_along(rho_vals) |> set_names(rho_vals),
        ~ {
          OM <- .fit[[.x]][["wi"]]
          dimnames(OM) <- dimnames(.S_train)
          OM
        })
  },
  CLIME = function(.fit, .S_train){
    
    map(seq_along(rho_vals) |> set_names(rho_vals),
        ~ {
          OM <- .fit[["icov"]][[.x]]
          if(all(OM == 0)){
            rho <- rho_vals[[.x]]
            OM <- diag(1/diag(.S_train))
          }
          
          dimnames(OM) <- dimnames(.S_train)
          OM
        })
  },
  
  SCIO = function(.fit, .S_train){
    
    map(seq_along(rho_vals) |> set_names(rho_vals),
        ~ {
          OM <- .fit[[.x]][["w"]]
          dimnames(OM) <- dimnames(.S_train)
          OM
        })
  }
)


extract_S_train_estimate <- switch(
  params$algo,
  QUIC = function(.fit, .S_train){
    
    map(seq_along(rho_vals) |> set_names(rho_vals),
        ~ {
          S <- .fit[["W"]][,, .x ]
          dimnames(S) <- dimnames(.S_train)
          S
        })
  },
  glasso = function(.fit, .S_train){
    
    map(seq_along(rho_vals) |> set_names(rho_vals),
        ~ {
          S <- .fit[[.x]][["w"]]
          dimnames(S) <- dimnames(.S_train)
          S
        })
  },
  CLIME = function(.fit, .S_train){
    # need to compute inverse
    map(seq_along(rho_vals) |> set_names(rho_vals),
        ~ {
          OM <- .fit[["icov"]][[.x]]
          
          if(all(OM == 0)){
            S <- diag(diag(.S_train))
          } else{
            S <- solve(OM)
          }
          
          dimnames(S) <- dimnames(.S_train)
          S
        })
  },
  
  SCIO = function(.fit, .S_train){
    
    map(seq_along(rho_vals) |> set_names(rho_vals),
        ~ {
          OM <- .fit[[.x]][["w"]]
          
          if(any(diag(OM) == 0)){
            prob_diags <- which(diag(OM) == 0)
            diag(OM)[prob_diags] <- abs(rnorm(length(prob_diags), sd = sd(diag(OM))/10))
          }
          
          S <- solve(OM)
          dimnames(S) <- dimnames(.S_train)
          S
        })
  }
)


estimate_se <- function(.OM, .sf_valid){
  OM21 <- .OM[(nb_se + 1):(nb_se + nb_sf), 1:nb_se]
  OM11 <- .OM[1:nb_se, 1:nb_se]
  
  W <- - OM21 %*% solve(OM11) # based on the estimated precision matrix
  t(t(W) %*% t(.sf_valid$mat))
}


untransform_se_hat <- function(se_valid_hat_t, prev_transform){
  
  transform_rev(se_valid_hat_t,
                prev_transform$parameters)
}

extract_adj_mat <- function(OM){
  OM[startsWith(rownames(OM), "SE_"),
     !startsWith(colnames(OM), "SE_")]
}


# get_coefs_from_OM <- function(OM){
#   
#   OM[startsWith(rownames(OM), "SE_"),
#      !startsWith(colnames(OM), "SE_")] |>
#     as.data.frame() |>
#     rownames_to_column("event_id_percount") |>
#     pivot_longer(cols = -event_id_percount,
#                  names_to = "sf_id",
#                  values_to = "coefficient") |>
#     mutate(event_id = str_extract(event_id_percount, "^SE_[0-9]+")) |>
#     select(-event_id_percount) |>
#     group_by(sf_id, event_id) |>
#     summarize(coefficient = max(coefficient),
#               .groups = "drop") |>
#     left_join(all_interactions_by_event,
#               by = c("event_id", "sf_id" = "sf_tx")) |>
#     mutate(literature = replace_na(literature, FALSE))
# }












