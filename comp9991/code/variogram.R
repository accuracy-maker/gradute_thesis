# Variogram Function
# Autor: Haitao Gao
# Date: 11st July 2025

library(parallel)


cal_vario_two_params <- function(coord, params, ncores = NULL) {
  lambda = params[1]
  vartheta = params[2]
  
  # assert the coordinates is a matrix
  if (!is.matrix(coord)) {
    coord = matrix(coord, nrow = 1) # convert it to 1-row matrix, where the number of columns equals the length of the vector
  }
  
  nrows = nrow(coord)
  
  # Edge case: single location
  if (nrows == 1) {
    val = (sqrt(sum(coord[1,]^2)) / lambda)^vartheta # norm2 coord: sqrt(x^2 + y^2 + z^2)
    return(val)
  }
  
  # define variogram between two location
  variogram_btw_two_locs <- function(coord) {
    # if the input is a vector, calculate the distance respect to origin
    if (is.vector(coord)) {
      val <- (sqrt(sum(coord^2)) / lambda)^vartheta
    }
    # else calculate the distance between two locations
    else {
      val <- (sqrt(sum((coord[1,] - coord[2,])^2)) / lambda)^vartheta
    }
    return(val)
  }
  
  # compute all pairs of index: all combinations (i,j), i < j, without duplication
  all_distrance_pair <- Rfast::comb_n(1:nrows, 2) # 2-row matrix and each column represents one pair. e.g. (1,2), (1,3), (1,4), (2,3)..
  
  # calculate two vectors:
  # 1. r(h_ij)
  # 2. r(0)
  
  if (is.null(ncores)) {
    gamma_vec <- unlist(lapply(1:ncol(all_distrance_pair), function(i) {
      idx = all_distrance_pair[, i]; variogram_btw_two_locs(coord[idx, ])
    }))
    gamma_origin = unlist(lapply(1:nrows, function(i) variogram_btw_two_locs(coord[i, ])))
  }
  else {
    gamma_vec = unlist(parallel::mclapply(1:ncol(all_distrance_pair), function(i) {
      idx = all_distrance_pair[, i]; variogram_btw_two_locs(coord[idx, ])
    }, mc.cores = ncores, mc.preschedule = TRUE))
    gamma_origin = unlist(parallel::mclapply(1:nrows, function(i) variogram_btw_two_locs(coord[i, ]), mc.cores = ncores, mc.preschedule = TRUE))
  }
  
  # cov(i,j) = 2r(0) - r(h_ij)
  cov_mat <- matrix(0, nrows, nrows) # create a [n, n] matrix
  cov_mat[t(all_distrance_pair)] <- gamma_vec # t:transpose; cov(i,j) = r(h_ij)
  cov_mat[t(all_distrance_pair[2:1, ])] <- gamma_vec # fill the symmetry part
  # Now, the cov(i,j) = cov(j,i) = r(h_ij)
  
  # Add the r(0_i) and r(0_j) to each element
  cov_mat <- matrix(gamma_origin, nrows, nrows, byrow = TRUE) +
    matrix(gamma_origin, nrows, nrows, byrow = FALSE) -
    cov_mat
  
  return(cov_mat + .Machine$double.eps * diag(nrows))
}


cal_vario_four_params <- function(coord, params, ncores=NULL){
  
  # define parameters
  lambda <- params[1]
  vartheta <- params[2]
  zeta <- params[3]
  m <- params[4]
  
  if(!is.matrix(coord)){coord <- matrix(coord,nrow=1)}
  
  # define the anisotropy matrix
  # cos(zeta) -sin(zeta)
  # m*sin(zeta) m*cos(zeta)
  
  V <- matrix(c(cos(zeta),m*sin(zeta),-sin(zeta),m*cos(zeta)), nrow = 2, ncol = 2)
  
  nrows <- nrow(coord)
  
  Omega = t(V) %*% V # Omega = V^T * V
  
  # Edge Case: single location
  if(nrows == 1){
    val <- (sqrt(c(coord[1,] %*% Omega %*% coord[1,])) / lambda)^vartheta
    return (val)
  }
  
  # vario function for two locations
  variogram_btw_two_locs_four_params <- function(coord){
    # if coord is not a matrix, compute r(0_i)
    if(!is.matrix(coord)){
      val <- (sqrt(c(coord %*% Omega %*% coord)) / lambda)^vartheta
    }
    # compute the pair vario r(h_ij)
    else{
      h = coord[1,] - coord[2,]
      val <- (sqrt(c(h %*% Omega %*% h)) / lambda)^vartheta
    }
    return (val)
  }


  # store all pair index
  all_pair_idx <- Rfast::comb_n(1:nrows, 2)
  
  # compute r(h)
  if(is.null(ncores)){
    # r(h_ij) when i != j
    gamma_vec <- unlist(lapply(1:ncol(all_pair_idx), function(i) {idx = all_pair_idx[,i]; variogram_btw_two_locs_four_params(coord[idx,])}))
    
    # r(0)
    gamma_origin <- unlist(lapply(1:nrows, function(i) variogram_btw_two_locs_four_params(coord[i,])))
  }
  else{
    gamma_vec <- unlist(mclapply(1:ncol(all_pair_idx), function(i) {idx = all_pair_idx[,i]; variogram_btw_two_locs_four_params(coord[idx,])}, mc.cores = ncores, mc.preschedule = TRUE))
    gamma_origin <- unlist(mclapply(1:nrows, function(i) variogram_btw_two_locs_four_params(coord[i,]), mc.cores = ncores, mc.preschedule = TRUE))
  }
  
  # compute the covariance
  cov_mat = matrix(0,nrows,nrows)
  cov_mat[t(all_pair_idx)] <- gamma_vec
  cov_mat[t(all_pair_idx[2:1,])] <- gamma_vec
  cov_mat = matrix(gamma_origin, nrows, nrows, byrow=TRUE) + matrix(gamma_origin, nrows, nrows, byrow=FALSE) - cov_mat
  
  return(cov_mat + .Machine$double.eps * diag(nrows))
}
  

###### Update lambda only
cal_vario_one_param <- function(coord, lambda, ncores = NULL) {
  
  # fixed parameters
  vartheta <- 1
  zeta <- -0
  m <- 1
  
  if (!is.matrix(coord)) {
    coord <- matrix(coord, nrow = 1)
  }
  
  # define the anisotropy matrix
  V <- matrix(c(cos(zeta), m * sin(zeta), -sin(zeta), m * cos(zeta)), nrow = 2, ncol = 2)
  
  nrows <- nrow(coord)
  Omega <- t(V) %*% V  # Omega = V^T * V
  
  # Edge Case: single location
  if (nrows == 1) {
    val <- (sqrt(c(coord[1, ] %*% Omega %*% coord[1, ])) / lambda)^vartheta
    return(val)
  }
  
  # vario function for two locations
  variogram_btw_two_locs_fixed <- function(coord) {
    if (!is.matrix(coord)) {
      val <- (sqrt(c(coord %*% Omega %*% coord)) / lambda)^vartheta
    } else {
      h <- coord[1, ] - coord[2, ]
      val <- (sqrt(c(h %*% Omega %*% h)) / lambda)^vartheta
    }
    return(val)
  }
  
  # store all pair index
  all_pair_idx <- Rfast::comb_n(1:nrows, 2)
  
  # compute r(h)
  if (is.null(ncores)) {
    gamma_vec <- unlist(lapply(1:ncol(all_pair_idx), function(i) {
      idx <- all_pair_idx[, i]
      variogram_btw_two_locs_fixed(coord[idx, ])
    }))
    
    gamma_origin <- unlist(lapply(1:nrows, function(i) variogram_btw_two_locs_fixed(coord[i, ])))
  } else {
    gamma_vec <- unlist(parallel::mclapply(1:ncol(all_pair_idx), function(i) {
      idx <- all_pair_idx[, i]
      variogram_btw_two_locs_fixed(coord[idx, ])
    }, mc.cores = ncores, mc.preschedule = TRUE))
    
    gamma_origin <- unlist(parallel::mclapply(1:nrows, function(i) {
      variogram_btw_two_locs_fixed(coord[i, ])
    }, mc.cores = ncores, mc.preschedule = TRUE))
  }
  
  # compute the covariance
  cov_mat <- matrix(0, nrows, nrows)
  cov_mat[t(all_pair_idx)] <- gamma_vec
  cov_mat[t(all_pair_idx[2:1, ])] <- gamma_vec
  cov_mat <- matrix(gamma_origin, nrows, nrows, byrow = TRUE) +
    matrix(gamma_origin, nrows, nrows, byrow = FALSE) - cov_mat
  
  return(cov_mat + .Machine$double.eps * diag(nrows))
}



########## Update lambda and vartheta
cal_vario_lv <- function(coord, lambda, vartheta, ncores = NULL) {
  
  # fixed parameters
  zeta <- -0.73
  m <- 0.96
  
  if (!is.matrix(coord)) {
    coord <- matrix(coord, nrow = 1)
  }
  
  # define the anisotropy matrix
  V <- matrix(c(cos(zeta), m * sin(zeta), -sin(zeta), m * cos(zeta)), nrow = 2, ncol = 2)
  nrows <- nrow(coord)
  Omega <- t(V) %*% V  # Omega = V^T * V
  
  # Edge Case: single location
  if (nrows == 1) {
    val <- (sqrt(c(coord[1, ] %*% Omega %*% coord[1, ])) / lambda)^vartheta
    return(val)
  }
  
  # variogram function between two locations
  variogram_btw_two_locs_fixed <- function(coord) {
    if (!is.matrix(coord)) {
      val <- (sqrt(c(coord %*% Omega %*% coord)) / lambda)^vartheta
    } else {
      h <- coord[1, ] - coord[2, ]
      val <- (sqrt(c(h %*% Omega %*% h)) / lambda)^vartheta
    }
    return(val)
  }
  
  # all pairs of indices
  all_pair_idx <- Rfast::comb_n(1:nrows, 2)
  
  # compute r(h)
  if (is.null(ncores)) {
    gamma_vec <- unlist(lapply(1:ncol(all_pair_idx), function(i) {
      idx <- all_pair_idx[, i]
      variogram_btw_two_locs_fixed(coord[idx, ])
    }))
    
    gamma_origin <- unlist(lapply(1:nrows, function(i) variogram_btw_two_locs_fixed(coord[i, ])))
  } else {
    gamma_vec <- unlist(parallel::mclapply(1:ncol(all_pair_idx), function(i) {
      idx <- all_pair_idx[, i]
      variogram_btw_two_locs_fixed(coord[idx, ])
    }, mc.cores = ncores, mc.preschedule = TRUE))
    
    gamma_origin <- unlist(parallel::mclapply(1:nrows, function(i) {
      variogram_btw_two_locs_fixed(coord[i, ])
    }, mc.cores = ncores, mc.preschedule = TRUE))
  }
  
  # construct covariance matrix
  cov_mat <- matrix(0, nrows, nrows)
  cov_mat[t(all_pair_idx)] <- gamma_vec
  cov_mat[t(all_pair_idx[2:1, ])] <- gamma_vec
  cov_mat <- matrix(gamma_origin, nrows, nrows, byrow = TRUE) +
    matrix(gamma_origin, nrows, nrows, byrow = FALSE) - cov_mat
  
  return(cov_mat + .Machine$double.eps * diag(nrows))
}


