transform_npn_shrinkage <- function(mat, parameters = NULL){
  huge::huge.npn(mat, verbose = FALSE)
}


# Takes a matrix (of SF or PSI in training set), transform it and
# also return the transformation parameters
transform_zscore <- function(mat, parameters = NULL){
  if(!is.null(parameters)){
    stopifnot(identical(names(parameters),
                        c("means", "sds")))
    params_given <- TRUE
  } else{
    parameters <- list()
    params_given <- FALSE
  }
  
  if(!params_given) parameters$means <- colMeans(mat, na.rm = TRUE)
  mat <- sweep(mat, 2L, parameters$means)
  
  
  if(!params_given) parameters$sds <- matrixStats::colSds(mat, na.rm = TRUE)
  mat <- sweep(mat, 2L, parameters$sds, `/`)
  list(mat = mat, parameters = parameters)
}

reverse_transform_zscore <- function(mat_trans, parameters = NULL){
  stopifnot(identical(names(parameters),
                      c("means", "sds")))
  
  mat <- sweep(mat_trans, 2L, parameters$sds, "*")
  mat <- sweep(mat, 2L, parameters$means, "+")
  
  col_all_NaN <- apply(mat, 2L, \(col) all(is.nan(col)))
  mat[,col_all_NaN] <- 0
  
  mat
}

##  Tests ----

# # test that Zscore does same as scale()
# all.equal(transform_zscore(mat_psi_train)$mat,
#           scale(mat_psi_train, center = TRUE, scale = TRUE),
#           check.attributes = FALSE)
# 
# # test that reverse reverses
# all.equal(mat_psi,
#           mat_psi |> transform_zscore() |> do.call(reverse_transform_zscore, args = _))
# 
# # test that we can reuse parameters
# bb <- transform_zscore(mat_psi_train)
# aa <- transform_zscore(mat_psi_train, bb$parameters)
# all.equal(bb$mat, aa$mat)
# 
# # test that we can transform a single row
# bb <- transform_zscore(mat_psi_train)
# aa <- transform_zscore(mat_psi_test[5,,drop=FALSE], bb$parameters)$mat
# aa
# 
# # test that our parameters are taken into account
# rndm_mat_tests <- matrix(1:3, nrow = 1)
# all.equal(transform_zscore(rndm_mat_tests,
#                            parameters = list(means = 3:1,
#                                              sds = rep(1,3)))$mat,
#           matrix(c(1-3, 2-2, 3-1), nrow = 1))
# 
# all.equal(transform_zscore(rndm_mat_tests,
#                            parameters = list(means = rep(1,3),
#                                              sds = 1:3))$mat,
#           matrix(c((1-1)/1, (2-1)/2, (3-1)/3), nrow = 1))



