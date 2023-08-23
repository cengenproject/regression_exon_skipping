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


loss_tong_cv_I <- function(Y, OM){
  D <- diag(diag(OM))
  
  norm(solve(D) %*% OM %*% t(Y),
             type = "F")^2
}
