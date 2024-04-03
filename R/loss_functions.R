logdet <- function(x){
  dt <- determinant(x, logarithm = TRUE)
  if(dt$sign != 1) return(-Inf)
  dt$modulus |> as.numeric()
}

loss_loglik <- function(Sts, OMtr){
  stopifnot(all(nrow(Sts) == nrow(OMtr),
                nrow(Sts) == ncol(OMtr),
                ncol(Sts) == nrow(Sts)))
  
  loglik <- -sum(Sts * OMtr) + logdet(OMtr) - nrow(Sts)*log(2*pi)
  loglik/2
}

loss_kl <- function(Sts, OMtr){
  stopifnot(all(nrow(Sts) == nrow(OMtr),
                nrow(Sts) == ncol(OMtr),
                ncol(Sts) == nrow(Sts)))
  
  0.5 * (sum(Sts * OMtr) - logdet(Sts %*% OMtr) - nrow(Sts))
}

# Frobenius loss. Note when na.rm=TRUE, can only compare result for same set of NA
loss_frob <- function(Sts, Str, na.rm = TRUE){
  stopifnot(all(nrow(Sts) == nrow(Str),
                nrow(Sts) == ncol(Str),
                ncol(Sts) == nrow(Sts)))
  
  # norm(Sts - Str, type = "F") # equivalent when no NA
  sqrt(sum( (Sts - Str)^2, na.rm = na.rm))
}

# Test Frobenius norm different methods
# A <- matrix(rnorm(12),nrow = 4)
# 
# A[3,2] <- NA
# A
# 
# norm(A, type = "F")
# 
# sqrt(sum(diag(t(A) %*% A)))
# 
# N <- 0
# for(i in seq_len(nrow(A))){
#   for(j in seq_len(ncol(A))){
#     N <- N + A[i,j]^2
#   }
# }
# sqrt(N)
# 
# sqrt(sum(A^2, na.rm = TRUE))



# matrix multiplication with na.rm = TRUE
# roughly equal to %*% when no NAs
# `%*na%` <- function(A, B){
#   stopifnot(ncol(A) == nrow(B))
#   
#   C <- matrix(nrow = nrow(A),
#               ncol = ncol(B))
#   
#   for(i in seq_len(nrow(C))){
#     for(j in seq_len(ncol(C))){
#       
#       C[i,j] <- sum(A[i,] * B[,j], na.rm = TRUE)
#       
#     }
#   }
#   C
# }

# more efficient version
`%*na%` <- function(A, B){
  A[is.na(A)] <- 0
  B[is.na(B)] <- 0
  A %*% B
}


loss_quad <- function(Sts, OMtr){
  stopifnot(all(nrow(Sts) == nrow(OMtr),
                nrow(Sts) == ncol(OMtr),
                ncol(Sts) == nrow(Sts)))
  
  sum(diag( ( Sts %*na% OMtr - diag(nrow(Sts)) )^2 ))
}



loss_tong_cv_I <- possibly(function(Y, OM){
  D <- diag(diag(OM))
  
  norm(solve(D) %*% OM %*% t(Y),
       type = "F")^2
},
otherwise = Inf)





perdiag <- function(nrow){
  K <- matrix(0L, nrow = nrow, ncol = nrow)
  for(i in 1:nrow){
    K[i,(nrow-i+1)] <- 1L
  }
  K
}

modif_chol <- function(M){
  stopifnot((N <- nrow(M)) == ncol(M))
  K <- perdiag(N)
  
  C <-  t(chol(K %*% M %*% K))
  # stopifnot(all.equal(K %*% M %*% K,
  #                     C %*% t(C)))
  
  S <- diag(diag(C))
  
  L <- C %*% solve(S)
  
  # Additional checks
  # D0 <- S^2
  # 
  # stopifnot(all.equal(K %*% M %*% K,
  #                     L %*% D0 %*% t(L)))
  
  
  # get to the same form as paper
  T <- t(K %*% L %*% K)
  
  T
  
  # Can be useful for checks, not used here:
  # D <- matrix(0,nrow(D0),ncol(D0)); diag(D) <- 1/rev(diag(D0))
  # 
  # stopifnot(all.equal(M,
  #                     t(T) %*% solve(D) %*% T,
  #                     check.attributes = FALSE))
  # list(T = T, D = D)
}


loss_tong_cv_II <- possibly(
  function(Y, OM){
    mat_T <- modif_chol(OM)
    
    norm(mat_T %*% t(Y), type = "F")^2
  },
  otherwise = Inf)


# Proportion of entries that are 0
mat_sparsity <- function(mat){
  
  sum(mat == 0, na.rm = TRUE)/length(mat)
}


# Fitting a power law to the connectivity, return Rsquare
mat_power_law <- function(adj){

  connectivity <- colSums(adj != 0)
  
  tab <- table(connectivity)
  if(length(tab) == 1L) return(NA_real_)
  
  k <- as.numeric(names(tab))[-1]
  p_k <- as.numeric(tab)[-1]
  
  cor(log10(k), log10(p_k))^2
}


compute_metric_rsquared <- function(measured, estimated){
  
  lm(as.numeric(estimated) ~ as.numeric(measured)) |>
    summary() |>
    (\(x) x[["adj.r.squared"]])()
}





frac_explained_var <- function(resid, measured, na.rm = FALSE){
  SSerr <- rowSums(resid^2, na.rm = na.rm)
  SStot <- apply(measured, 1, \(.x) sum((.x - mean(.x, na.rm = na.rm))^2, na.rm = na.rm))
  
  pmax(1 - SSerr/SStot, 0)
}


compute_metric_FEV <- function(measured, estimated){
  mat_residuals <- estimated - measured
  sum_abs_residuals <- sum(abs(mat_residuals), na.rm = TRUE)
  
  FEV <- frac_explained_var(mat_residuals, measured, na.rm = TRUE)
  mean_FEV <- mean(FEV)
  
  mean_FEV
}




impute_median <- function(mat){
  
  apply(mat, 2L, \(x) {x[is.na(x)] <- median(x, na.rm = TRUE); x})
}

impute_knn <- function(mat){
  old_seed <- .Random.seed
  res <- impute::impute.knn(mat, k = params$knn_k)[["data"]]
  .Random.seed <- old_seed
  res
}






compute_TPR <- function(adj){
  
  rows <- intersect(rownames(adj), rownames(mat_interactions_lit))
  cols <- intersect(colnames(adj), colnames(mat_interactions_lit))
  
  nb_TP <- sum( (adj[rows, cols] != 0) * mat_interactions_lit[rows, cols] )
  nb_T <- sum(mat_interactions_lit[rows,cols])
  
  nb_TP / nb_T
}

compute_FPR <- function(adj){
  
  rows <- intersect(rownames(adj), rownames(mat_interactions_lit))
  cols <- intersect(colnames(adj), colnames(mat_interactions_lit))
  
  nb_FP <- sum( (adj[rows, cols] != 0) * (!mat_interactions_lit[rows, cols]) )
  nb_F <- sum( ! mat_interactions_lit[rows,cols] )
  
  nb_FP / nb_F
}



reconstruct_psi_from_counts <- function(mat_N){
  stopifnot(identical(str_split_i(colnames(mat_N)[c(TRUE,FALSE)], "\\.", 1),
                      str_split_i(colnames(mat_N)[c(FALSE,TRUE)], "\\.", 1)))
  
  stopifnot(identical(str_replace(colnames(mat_N)[c(TRUE,FALSE)],
                                  "\\.Nincl", "\\.Nexcl"),
                      colnames(mat_N)[c(FALSE,TRUE)]))
  
  mat_Nincl <- mat_N[,c(TRUE,FALSE)]
  colnames(mat_Nincl) <- str_remove(colnames(mat_Nincl), "\\.Nincl$")
  mat_Nexcl <- mat_N[,c(FALSE,TRUE)]
  colnames(mat_Nexcl) <- str_remove(colnames(mat_Nexcl), "\\.Nexcl$")
  
  stopifnot(all.equal(dimnames(mat_Nincl), dimnames(mat_Nexcl)))
  
  
  mat_Nincl/(mat_Nincl + mat_Nexcl)
}



reconstruct_adj_psi <- function(adj_cnt){
  adj_incl <- adj_cnt[c(TRUE, FALSE), ]
  rownames(adj_incl) <- str_remove(rownames(adj_incl), "\\.Nincl$")
  
  adj_excl <- adj_cnt[c(FALSE, TRUE), ]
  rownames(adj_excl) <- str_remove(rownames(adj_excl), "\\.Nexcl$")
  
  stopifnot(all.equal(dimnames(adj_incl),
                      dimnames(adj_excl)))
  
  adj_psi <- adj_incl - adj_excl
  adj_psi
}

# average sparsity by event
adj_sparsity <- function(adj_psi){
  
  rowMeans(adj_psi != 0) |> mean()
}









