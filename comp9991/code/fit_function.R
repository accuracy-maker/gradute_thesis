# Fit Function
# Author: Haitao Gao
# Date: 13th July 2025

library(parallel)
################## Intensity Function ##############
# the expression for kappa(x), which is the intensity function of log-skewed Brown-Resnick model
# W(s) = exp{Y(s) - a(s)}, where Y(s) is a centred skew-normal process with scale matrix \Sigma and slant paramter \alpha.
# PARAMS:
# 1. x: the variable
# 2. params: parameters: params[[1]] = \Sigma, params[[2]] = \alpha
# 3. alpha_para: the parameters about alpha
# 4. log: logrithm or not

intensity_logskew_constraint <- function(x, params, is_alpha=TRUE, log=TRUE){
  oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
  set.seed(747380)
  
  # extract Sigma and alpha
  Sigma <- params[[1]] # this is computed by variogram function in "variogram.R" Script
  
  if(!is.matrix(x)) x <- matrix(x, nrow = 1) # convert x to matrix if it's a vector
  n = ncol(x)
  
  # Edge case: 1D data (single location)
  if(n == 1) return (1 / (x^2))
  
  # pre-compute the core matrix terms used for Gaussian terms
  # omega^2 = Var(X_i) = diag(Sigma)
  # cholesky decomposition
  # inverse of Sigma
  # log(Sigma)
  
  omega_pow_2 <- diag(Sigma)
  chol_Sigma <- chol(Sigma)
  inv_Sigma <- chol2inv(chol_Sigma)
  logdet_Sigma <- sum(log(diag(chol_Sigma)))*2
  
  # if we have alpha part (skewed) otherwise it's the standard Brown-Resnick model
  if(is_alpha){
    alpha <- params[[2]]
    # derive delta by alpha
    delta <- c(Sigma %*% alpha) / sqrt(c(1+alpha %*% Sigma %*% alpha))
  }
  # given delta directly
  else{
    delta <- params[[2]]
    # derive alpha by delta
    alpha = c(1 - delta %*% inv_Sigma %*% delta)^(-1/2) * c(inv_Sigma %*% delta)
  }
  
  # pre-compute the constant
  # constant a(s)
  # constant q
  
  a = log(2) + pnorm(delta, log.p = TRUE)
  q = rowSums(inv_Sigma)
  sum_q = sum(q)
  
  q_mat <- matrix(q, n, n, byrow = TRUE)
  
  # transform data and defnie skew-normal part
  # x_i^o = logx_i + a_i
  x_log <- log(x)
  x_circ <- x_log + matrix(a, nrow = nrow(x), ncol = n, byrow = TRUE)
  
  # compute tilt parameter
  # tau = alpha^T * (x_circ + oemga^2 / 2)
  tau_tilde <- apply(x_circ, 1, function(x_i) sum(alpha * (x_i + omega_pow_2)))
  
  # Cnetering matrix
  # A projection matrix: project into orthogonal subspace
  # A = inv_Sigma - q * q^T / sum(q)
  A = inv_Sigma - q %*% t(q) / sum_q
  
  # full log-density expression
  val = -(n-3)/2 * log(2) - (n-1)/2*log(pi)-1/2*logdet_Sigma - 1/2*log(sum_q) - 1/2 * (sum(q*omega_pow_two)-1)/sum_q - 1/8*c(omega_pow_2 %*% A %*% omega_pow_2) 
  val = val - rowSums(x_log) - 1/2 * apply(x_circ,1,function(x_i) c(x_i %*% A %*% x_i) + sum(x_i * (2*q/sum_q + c(A %*% omega_pow_2)))) + pnorm(tau_tilde,log.p=TRUE)
  assign(".Random.seed", oldSeed, envir=globalenv())
  
  if(log){
    return (val)
  }
  else{
    return(exp(val))
  }
}

############# Alpha Function ##############
alpha_fun <- function(params, b_mat=basis){
  alpha <- c(params %*% t(b_mat))
  return (alpha)
}


################# Function ###############
# objective: find the parameters by optimise the neg-likelihood function
# INPUT: 
# 1. data: extreme values above one specific threshold u
# 2. coords: coordinates of locations used to capture the spatial dependency
# 3. init: initial values of parameters
# 4. maxiter: maximum iterations
# 5. FUN: which variogram function we use: two parameters or four parameters
# 6. ncore: how many CPU cores to use
# 7. optimiser: "Nelder-Mead"
# 8. lb: lower bound
# 9. ub: upper bound
# 10: 


fit_model <- function(data, coords, init, params_fixed = NULL,
                      maxit = 100, vario_func = NULL, basis = NULL,
                      alpha_func = NULL, ncores = NULL, partial = TRUE,
                      optimiser = "Nelder-Mead", lb = NULL, ub = NULL,
                      hessian = FALSE, opt = TRUE, step2 = TRUE,
                      idx_para = 1:2) {
  
  t0 <- proc.time()
  
  if (is.null(lb)) lb <- rep(-Inf, length(init))
  if (is.null(ub)) ub <- rep(Inf, length(init))
  if (is.null(params_fixed)) params_fixed <- rep(FALSE, length(init))
  fixed <- params_fixed
  fixed2 <- fixed
  
  object_func <- function(params, opt = TRUE, ncores = ncores) {
    if (any(params < lb[!fixed2]) | any(params > ub[!fixed2])) {return(Inf)}
    
    params2 <- init
    params2[!fixed2] <- params
    
    param_sigma <- params2[idx_para]
    param_alpha <- params2[-idx_para]
    
    sigma <- vario_func(coords, param_sigma, ncores)
    alpha <- alpha_func(param_alpha, basis)
    
    delta <- c(sigma %*% alpha) / sqrt(c(1 + alpha %*% sigma %*% alpha))
    omega2 <- diag(sigma)
    a <- log(2) + pnorm(delta, log.p = TRUE)
    
    computeScores <- function(i) {
      ind <- data[[i]][[1]]
      x <- data[[i]][[2]]
      n <- length(x)
      if (n == 1) return(1 / (x^2))
      sigma_i <- sigma[ind, ind]
      chol_sigma <- chol(sigma_i)
      inv_sigma <- chol2inv(chol_sigma)
      logdet_sigma <- sum(log(diag(chol_sigma))) * 2
      delta_i <- delta[ind]
      omega2_i <- omega2[ind]
      alpha_i <- c(1 - delta_i %*% inv_sigma %*% delta_i)^(-1/2) * c(inv_sigma %*% delta_i)
      a_i <- a[ind]
      q <- rowSums(inv_sigma)
      sum_q <- sum(q)
      sum_alpha <- sum(alpha)
      x_log <- log(x)
      x_circ <- x_log + a_i
      beta <- (1 + sum_alpha^2 / sum_q)^(-0.5)
      tau_tilde <- beta * sum((alpha_i - sum_alpha * q / sum_q) * (x_circ + omega2_i / 2)) +
        beta * sum_alpha / sum_q
      A <- inv_sigma - q %*% t(q) / sum_q
      val <- -(n - 3) / 2 * log(2) - (n - 1) / 2 * log(pi) - 0.5 * logdet_sigma -
        0.5 * log(sum_q) - 0.5 * (sum(q * omega2_i) - 1) / sum_q -
        0.125 * c(omega2_i %*% A %*% omega2_i)
      val <- val - x_log - 0.5 * (c(x_circ %*% A %*% x_circ) +
                                    sum(x_circ * (2 * q / sum_q + c(A %*% omega2_i)))) +
        pnorm(tau_tilde, log.p = TRUE)
      return(val)
    }
    
    if (!is.null(ncores)) {
      val <- unlist(parallel::mclapply(1:length(data), computeScores, mc.cores = ncores, mc.preschedule = TRUE))
    } else {
      val <- unlist(lapply(1:length(data), computeScores))
    }
    
    return(-mean(val))
  }
  
  if (opt) {
    opt_result <- optim(
      par = init[!fixed2],
      fn = object_func,
      method = optimiser,
      control = list(maxit = maxit, trace = FALSE, reltol = 1e-8),
      hessian = hessian,
      ncore = ncores
    )
    
    if (any(!fixed[-idx_para]) && step2 && !is.null(ncores)) {
      n_alpha <- sum(!fixed[-idx_para])
      ncores_2 <- ceiling(ncores / 4)
      
      if (n_alpha == 2) {
        a <- seq(0, 2 * pi, length.out = 4)
        a <- cbind(cos(a), sin(a))
      } else {
        a <- matrix(rnorm(ncores * n_alpha), ncol = n_alpha)
        a <- sweep(a, 1, sqrt(rowSums(a^2)), FUN = "/")
      }
      
      init[!fixed2] <- opt_result$par
      fixed2 <- fixed
      fixed2[idx_para] <- TRUE
      fixed2[-idx_para] <- fixed[-idx_para]
      init_list <- split(a, row(a))
      
      opt_result2 <- mcmapply(
        optim,
        par = init_list,
        MoreArgs = list(
          fn = object_func,
          method = optimiser,
          control = list(maxit = maxit, trace = FALSE, reltol = 1e-8),
          hessian = FALSE,
          ncore = ncores_2
        ),
        mc.cores = 4,
        mc.set.seed = FALSE,
        SIMPLIFY = FALSE
      )
      
      opt_values <- unlist(lapply(opt_result2, function(x) {
        tryCatch(x$value, error = function(e) Inf)
      }))
      opt_result <- opt_result2[[which.min(opt_values)]]
      
      init[!fixed2] <- opt_result$par
      fixed2 <- fixed
      
      opt_result <- optim(
        par = init[!fixed2],
        fn = object_func,
        method = optimiser,
        control = list(maxit = maxit, trace = FALSE, reltol = 1e-8),
        hessian = hessian,
        ncore = ncores
      )
      
      opt_result$others <- opt_result2
    }
    
  } else {
    opt_result <- list(par = init, value = object_func(init, opt = FALSE))
  }
  
  opt_result$time <- proc.time() - t0
  
  if (hessian) {
    opt_result$grad <- numDeriv::jacobian(object_func, opt_result$par, opt = FALSE)
    opt_result$hessian <- numDeriv::hessian(object_func, opt_result$par, opt = TRUE)
    opt_result$K <- var(opt_result$grad)
  }
  
  params2 <- init
  params2[!fixed2] <- opt_result$par
  opt_result$par <- params2
  
  return(opt_result)
}


################## Fit model by two risk functions ##############
fit_model_2 <- function(data, coords, init, params_fixed = NULL,
                            maxit = 100, vario_func = NULL, basis = NULL,
                            alpha_func = NULL, ncores = NULL,
                            optimiser = "L-BFGS-B", lb = NULL, ub = NULL,
                            hessian = FALSE, opt = TRUE, idx_para = 1:2) {
  
  t0 <- proc.time()
  if (is.null(lb)) lb <- rep(-Inf, length(init))
  if (is.null(ub)) ub <- rep(Inf, length(init))
  if (is.null(params_fixed)) params_fixed <- rep(FALSE, length(init))
  fixed2 <- params_fixed
  
  cat(">>> Optimizer:", optimiser, "\n")
  cat(">>> Using", ifelse(is.null(ncores), 1, ncores), "core(s)\n")
  
  # Split data by label
  data_0 <- Filter(function(x) x[[3]] == 0, data)
  data_1 <- Filter(function(x) x[[3]] == 1, data)
  data_2 <- Filter(function(x) x[[3]] == 2, data)
  
  object_func <- function(params, opt = TRUE, ncores = ncores) {
    if (any(params < lb[!fixed2]) | any(params > ub[!fixed2])) {return(Inf)}
    full_params <- init
    full_params[!fixed2] <- params
    
    n_total <- length(full_params) / 2
    set1 <- full_params[1:n_total]
    set2 <- full_params[(n_total + 1):(2 * n_total)]
    param_avg <- (set1 + set2) / 2
    
    compute_negloglik <- function(subdata, param_set) {
      param_sigma <- param_set[idx_para]
      param_alpha <- param_set[-idx_para]
      sigma <- vario_func(coords, param_sigma, ncores)
      alpha <- alpha_func(param_alpha, basis)
      
      delta <- c(sigma %*% alpha) / sqrt(c(1 + alpha %*% sigma %*% alpha))
      omega2 <- diag(sigma)
      a <- log(2) + pnorm(delta, log.p = TRUE)
      
      computeScores <- function(i) {
        ind <- subdata[[i]][[1]]
        x <- subdata[[i]][[2]]
        n <- length(x)
        if (n == 1) return(1 / (x^2))
        sigma_i <- sigma[ind, ind]
        chol_sigma <- chol(sigma_i)
        inv_sigma <- chol2inv(chol_sigma)
        logdet_sigma <- sum(log(diag(chol_sigma))) * 2
        delta_i <- delta[ind]
        omega2_i <- omega2[ind]
        alpha_i <- c(1 - delta_i %*% inv_sigma %*% delta_i)^(-1/2) * c(inv_sigma %*% delta_i)
        a_i <- a[ind]
        q <- rowSums(inv_sigma)
        sum_q <- sum(q)
        sum_alpha <- sum(alpha)
        x_log <- log(x)
        x_circ <- x_log + a_i
        beta <- (1 + sum_alpha^2 / sum_q)^(-0.5)
        tau_tilde <- beta * sum((alpha_i - sum_alpha * q / sum_q) * (x_circ + omega2_i / 2)) +
          beta * sum_alpha / sum_q
        A <- inv_sigma - q %*% t(q) / sum_q
        val <- -(n - 3) / 2 * log(2) - (n - 1) / 2 * log(pi) - 0.5 * logdet_sigma -
          0.5 * log(sum_q) - 0.5 * (sum(q * omega2_i) - 1) / sum_q -
          0.125 * c(omega2_i %*% A %*% omega2_i)
        val <- val - x_log - 0.5 * (c(x_circ %*% A %*% x_circ) +
                                      sum(x_circ * (2 * q / sum_q + c(A %*% omega2_i)))) +
          pnorm(tau_tilde, log.p = TRUE)
        return(val)
      }
      
      if (!is.null(ncores)) {
        val <- unlist(parallel::mclapply(1:length(subdata), computeScores, mc.cores = ncores, mc.preschedule = TRUE))
      } else {
        val <- unlist(lapply(1:length(subdata), computeScores))
      }
      return(sum(val))
    }
    
    # Evaluate NLL for three groups
    neglog_0 <- compute_negloglik(data_0, param_avg)
    neglog_1 <- compute_negloglik(data_1, set1)
    neglog_2 <- compute_negloglik(data_2, set2)
    
    return(-1 * mean(neglog_0 + neglog_1 + neglog_2))
  }
  
  if (opt) {
    free_init <- init[!fixed2]
    cat(">>> Initial free parameters:\n")
    print(free_init)
    
    if (optimiser == "L-BFGS-B") {
      opt_result <- optim(
        par = free_init,
        fn = object_func,
        method = "L-BFGS-B",
        lower = lb[!fixed2],
        upper = ub[!fixed2],
        control = list(maxit = maxit, trace = TRUE),
        hessian = hessian,
        ncore=ncores
      )
    } else {
      opt_result <- optim(
        par = free_init,
        fn = object_func,
        method = optimiser,
        control = list(maxit = maxit, trace = TRUE, reltol = 1e-8),
        hessian = hessian,
        ncore=ncores
      )
    }
  } else {
    opt_result <- list(par = init, value = object_func(init, opt = FALSE))
  }
  
  opt_result$time <- proc.time() - t0
  
  if (hessian) {
    opt_result$grad <- numDeriv::jacobian(object_func, opt_result$par, opt = FALSE)
    opt_result$hessian <- numDeriv::hessian(object_func, opt_result$par, opt = TRUE)
    opt_result$K <- var(opt_result$grad)
  }
  
  return(opt_result)
}






############### new version
# Updated function implementing Algorithm 1
# Renamed to fit_model_gaus for clarity
fit_model_gaus <- function(data, coords, init, 
                           maxit = 1000, 
                           vario_func = NULL, 
                           basis = NULL,
                           alpha_func = NULL,
                           ncores = NULL,
                           optimiser = "L-BFGS-B", 
                           lb = NULL, ub = NULL,
                           hessian = FALSE, opt = TRUE,
                           is_logskewed = FALSE,
                           mode = 0) {
  
  t0 <- proc.time()
  
  # Validate initial parameters based on mode
  if (mode == 0 && length(init) != 3) {
    stop("For mode 0, 'init' must be a vector of length 3 for c_1_inf, c_1, and c_inf.")
  } else if (mode %in% c(1, 2) && length(init) != 1) {
    stop("For mode 1 or 2, 'init' must be a vector of length 1.")
  }
  
  if (is.null(lb)) lb <- rep(-Inf, length(init))
  if (is.null(ub)) ub <- rep(Inf, length(init))
  
  cat(">>> Optimizer:", optimiser, "\n")
  cat(">>> Using", ifelse(is.null(ncores), 1, ncores), "core(s)\n")
  cat(">>> Training mode:", mode, "\n")
  
  # --- Data Subsetting based on Algorithm 1 ---
  # D_1 contains data with label = 0 or label = 1
  # D_inf contains data with label = 0 or label = 2
  data_for_S1 <- Filter(function(x) x[[3]] %in% c(0, 1), data)
  data_for_S_inf <- Filter(function(x) x[[3]] %in% c(0, 2), data)
  
  if (mode == 0) {
    cat(">>> Fitting c_1_inf, c_1, c_inf using D_1 and D_inf\n")
  } else if (mode == 1) {
    cat(">>> Fitting c_1 using D_1 (labels 0 & 1)\n")
  } else if (mode == 2) {
    cat(">>> Fitting c_inf using D_inf (labels 0 & 2)\n")
  } else {
    stop("Invalid mode. Must be 0, 1, or 2.")
  }
  
  # Objective function to be minimized (Total Negative Log-Likelihood)
  object_func <- function(params, opt = TRUE) {
    
    total_negloglik <- 0
    
    # --- Negative Log-Likelihood Calculation (reusable component) ---
    compute_negloglik <- function(subdata, lambda_val) {
      sigma <- vario_func(coords, lambda_val, ncores)
      alpha <- rep(0, nrow(basis)) # Assuming alpha is zero as per old code
      delta <- rep(0, nrow(sigma))
      omega2 <- diag(sigma)
      a <- log(2) + pnorm(delta, log.p = TRUE)
      
      computeScores <- function(i) {
        ind <- subdata[[i]][[1]]
        x <- subdata[[i]][[2]]
        n <- length(x)
        if (n == 1) return(1 / (x^2)) # Simplified case
        
        sigma_i <- sigma[ind, ind]
        chol_sigma <- chol(sigma_i)
        inv_sigma <- chol2inv(chol_sigma)
        logdet_sigma <- sum(log(diag(chol_sigma))) * 2
        
        delta_i <- delta[ind]
        omega2_i <- omega2[ind]
        a_i <- a[ind]
        
        q <- rowSums(inv_sigma)
        sum_q <- sum(q)
        x_log <- log(x)
        x_circ <- x_log + a_i
        
        beta <- 1
        tau_tilde <- beta * sum(x_circ + omega2_i / 2) / sqrt(sum_q)
        
        A <- inv_sigma - q %*% t(q) / sum_q
        val <- -(n - 3) / 2 * log(2) - (n - 1) / 2 * log(pi) - 0.5 * logdet_sigma -
          0.5 * log(sum_q) - 0.5 * (sum(q * omega2_i) - 1) / sum_q -
          0.125 * c(omega2_i %*% A %*% omega2_i)
        
        val <- val - x_log - 0.5 * (c(x_circ %*% A %*% x_circ) +
                                      sum(x_circ * (2 * q / sum_q + c(A %*% omega2_i)))) +
          pnorm(tau_tilde, log.p = TRUE)
        return(val)
      }
      
      if (!is.null(ncores)) {
        val <- unlist(parallel::mclapply(1:length(subdata), computeScores, mc.cores = ncores, mc.preschedule = TRUE))
      } else {
        val <- unlist(lapply(1:length(subdata), computeScores))
      }
      # The function computes log-likelihood, so we return the negative sum
      return(-sum(val))
    }
    
    # --- Calculate Total Negative Log-Likelihood based on Mode ---
    if (mode == 0) {
      # Params: c(c_1_inf, c_1, c_inf)
      c_1_inf <- params[1]
      c_1     <- params[2]
      c_inf   <- params[3]
      
      lambda_S1    <- c_1_inf + c_1
      lambda_S_inf <- c_1_inf + c_inf
      
      negloglik_D1 <- compute_negloglik(data_for_S1, lambda_S1)
      negloglik_D_inf <- compute_negloglik(data_for_S_inf, lambda_S_inf)
      
      total_negloglik <- negloglik_D1 + negloglik_D_inf
      
    } else if (mode == 1) {
      # Params: c(c_1)
      c_1 <- params[1]
      
      # Per algorithm, lambda for D1 is just c1
      lambda_S1 <- c_1
      total_negloglik <- compute_negloglik(data_for_S1, lambda_S1)
      
    } else if (mode == 2) {
      # Params: c(c_inf)
      c_inf <- params[1]
      
      # Per algorithm, lambda for D_inf is just c_inf
      lambda_S_inf <- c_inf
      total_negloglik <- compute_negloglik(data_for_S_inf, lambda_S_inf)
    }
    
    return(total_negloglik)
  }
  
  # Optimization
  if (opt) {
    cat(">>> Initial parameters:\n")
    print(init)
    
    opt_result <- optim(
      par = init,
      fn = object_func,
      method = optimiser,
      lower = lb,
      upper = ub,
      control = list(maxit = maxit, trace = TRUE),
      hessian = hessian
    )
  } else {
    opt_result <- list(par = init, value = object_func(init, opt = FALSE))
  }
  
  opt_result$time <- proc.time() - t0
  opt_result$mode <- mode
  
  # Add parameter names to the output for clarity
  if (mode == 0) {
    names(opt_result$par) <- c("c_1_inf", "c_1", "c_inf")
  } else if (mode == 1) {
    names(opt_result$par) <- "c_1"
  } else if (mode == 2) {
    names(opt_result$par) <- "c_inf"
  }
  
  if (hessian) {
    opt_result$grad <- numDeriv::jacobian(object_func, opt_result$par, opt = FALSE)
    opt_result$hessian <- numDeriv::hessian(object_func, opt_result$par, opt = TRUE)
    opt_result$K <- var(opt_result$grad)
  }
  
  return(opt_result)
}





