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
  T <- ncol(Y)
  N <- nrow(Y)
  
  # split
  N1 <- floor(.5*N)
  N2 <- N - N1
  
  Yp1 <- Y[1:N1,]
  Yp2 <- Y[N1 + 1:N2,]
  
  
  # partial precision matrix
  OM11_train <- OM[1:N1, 1:N1]
  OM21_train <- OM[(N1+1):(N1+N2), 1:N1]
  
  OM22_train <- OM[N1+(1:N2), N1+(1:N2)]
  OM12_train <- OM[1:N1, (N1+1):(N1+N2)]
  
  W_train1 <- - OM21_train %*% solve(OM11_train)
  W_train2 <- - OM12_train %*% solve(OM22_train)
  
  estimated_Yp1 <- t(W_train1) %*% Yp2
  estimated_Yp2 <- t(W_train2) %*% Yp1
  
  # Fit metrics
  Rsq1 <- cor(as.numeric(estimated_Yp1), as.numeric(Yp1))^2
  
  mae1 <- sum(abs(estimated_Yp1 - Yp1))/length(Yp1)
  
  sign_dpsi_matches1 <- sign(estimated_Yp1) == sign(Yp1)
  accuracy1 <- sum(sign_dpsi_matches1)/length(Yp1)
  
  
  
  Rsq2 <- cor(as.numeric(estimated_Yp2), as.numeric(Yp2))^2
  
  mae2 <- sum(abs(estimated_Yp2 - Yp2))/length(Yp2)
  
  sign_dpsi_matches2 <- sign(estimated_Yp2) == sign(Yp2)
  accuracy2 <- sum(sign_dpsi_matches2)/length(Yp2)
  
  tibble(Rsquare1 = Rsq1,
         mae1 = mae1,
         accuracy_sign1 = accuracy1,
         Rsquare2 = Rsq2,
         mae2 = mae2,
         accuracy_sign2 = accuracy2)
}
