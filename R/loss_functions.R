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

loss_frob <- function(Sts, Str){
  stopifnot(all(nrow(Sts) == nrow(Str),
                nrow(Sts) == ncol(Str),
                ncol(Sts) == nrow(Sts)))
  
  norm(Sts - Str, type = "F")
}

loss_quad <- function(Sts, OMtr){
  stopifnot(all(nrow(Sts) == nrow(OMtr),
                nrow(Sts) == ncol(OMtr),
                ncol(Sts) == nrow(Sts)))
  
  sum(diag( ( Sts %*% OMtr - diag(nrow(Sts)) )^2 ))
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


frac_explained_var <- function(resid, measured, na.rm = FALSE){
  SSerr <- rowSums(resid^2, na.rm = na.rm)
  SStot <- apply(measured, 1, \(.x) sum((.x - mean(.x, na.rm = na.rm))^2, na.rm = na.rm))
  
  pmax(1 - SSerr/SStot, 0)
}
