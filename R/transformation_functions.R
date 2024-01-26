

# NPN shrinkage ----


which_smallest_positive <- function(a){
  if(all(a <= 0 | is.na(a))){
    positive_min <- 0
  } else{
    positive_min <- min(a[a>=0], na.rm = TRUE)
  }
  
  which(a == positive_min)[[1]]
}

transform_one_value <- function(y1, x){
  if(is.na(y1)) return(length(x[!is.na(x)]) + 1)
  
  if(y1 < min(x, na.rm = TRUE)){
    right_neighbor <- which_smallest_positive( -(y1 - x))
    rank_x <- rank(x)
    
    return(rank_x[right_neighbor])
  } else if(y1 > max(x, na.rm = TRUE)){
    left_neighbor <- which_smallest_positive(y1 - x)
    rank_x <- rank(x)
    
    return(rank_x[left_neighbor])
  }
  
  left_neighbor <- which_smallest_positive(y1 - x)
  right_neighbor <- which_smallest_positive( -(y1 - x))
  
  rank_x <- rank(x)
  
  if(x[right_neighbor] == x[left_neighbor]) return(rank_x[left_neighbor])
  
  rank_x[left_neighbor] + (y1 - x[left_neighbor])*(rank_x[right_neighbor] - rank_x[left_neighbor])/(x[right_neighbor] - x[left_neighbor])
}

project_rank <- function(vec_to_transform, original_vec){
  map_dbl(vec_to_transform, transform_one_value, original_vec)
}



transform_npn_shrinkage <- function(mat, parameters = NULL){
  
  if(!is.null(parameters)){
    stopifnot(identical( names(parameters),
                         c("reference_mat", "sd_first_col")
                        ))
    stopifnot(identical( colnames(mat),
                         colnames(parameters$reference_mat)
                        ))
    params_given <- TRUE
  } else{
    parameters <- list()
    params_given <- FALSE
  }
  
  if(! params_given){
    mat_trans <- apply(mat, 2, rank)
    mat_trans <- mat_trans/(nrow(mat_trans) + 1)
  } else{
    mat_trans <- matrix(nrow = nrow(mat),
                        ncol = ncol(mat))
    for(col in seq_len(ncol(mat))){
      mat_trans[,col] <- project_rank(mat[,col], parameters$reference_mat[,col])
    }
    dimnames(mat_trans) <- dimnames(mat)
    mat_trans <- mat_trans/(nrow(parameters$reference_mat) + 1)
  }
  
  mat_trans <- qnorm(mat_trans)
  
  if(params_given){
    sd_first_col <- parameters$sd_first_col
  } else{
    sd_first_col <- sd(mat_trans[, 1])
  }
  mat_trans = mat_trans/sd_first_col
   
  list(mat = mat_trans, parameters = list(reference_mat = mat, sd_first_col = sd_first_col))
}





## Tests ----

# # test that npn does same as {huge}
# all.equal(transform_npn_shrinkage(mat_psi)$mat,
#           huge::huge.npn(mat_psi, verbose = FALSE))
# 
# all.equal(transform_npn_shrinkage(mat_sf)$mat,
#           huge::huge.npn(mat_sf, verbose = FALSE))
# 
# 
# 
# # test that we can reuse parameters
# bb <- transform_npn_shrinkage(mat_sf_train)
# aa <- transform_npn_shrinkage(mat_sf_train, bb$parameters)
# all.equal(bb$mat, aa$mat)
# # with na: not an exact match, the huge function gives incresing ranks to each NA, I give them all the same rank
# bb <- transform_npn_shrinkage(mat_psi_train)
# aa <- transform_npn_shrinkage(mat_psi_train, bb$parameters)
# all.equal(bb$mat, aa$mat)
# plot(aa$mat,bb$mat, col = c("black", "darkred")[is.na(mat_psi_train) +1])
# 
# 
# 
# # test that we can transform a single row
# bb <- transform_npn_shrinkage(mat_psi_train)
# aa <- transform_npn_shrinkage(mat_psi_test[5,,drop=FALSE], bb$parameters)$mat
# aa[,1:3]
# 
# plot(log10(mat_psi_train), bb$mat)
# points(log10(mat_psi_test[5,]), aa, col = 'red')
# 
# 
# # test that our parameters are taken into account
# rndm_mat_tests <- matrix(1:3, nrow = 1)
# all.equal(transform_npn_shrinkage(rndm_mat_tests,
#                                   parameters = list(reference_mat = matrix(rep(1:3, 3), nrow = 3),
#                                                     sd_first_col = 1))$mat,
#           matrix(c(qnorm(1/4), qnorm(2/4), qnorm(3/4)),
#                  nrow = 1))
# 
# all.equal(transform_npn_shrinkage(rndm_mat_tests,
#                                   parameters = list(reference_mat = matrix(rep(1:3, 3), nrow = 3),
#                                                     sd_first_col = 2.3))$mat,
#           matrix(c(qnorm(1/4), qnorm(2/4), qnorm(3/4))/2.3,
#                  nrow = 1))
# 
# 
# all.equal(transform_npn_shrinkage(rndm_mat_tests,
#                                   parameters = list(reference_mat = matrix(rep(c(1,1,3), 3), nrow = 3),
#                                                     sd_first_col = 1))$mat,
#           matrix(c(qnorm(1.5/4), qnorm((1.5+(3-1.5)/2)/4), qnorm(3/4)),
#                  nrow = 1))






# Z-score ----


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



