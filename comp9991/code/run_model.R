# Load dependencies and custom functions
source("load_data.R")
source("variogram.R")
source("fit_function.R")

library(pbapply)
library(parallel)
library(mvtnorm)
library(evd)
library(Rfast)

# Set seed for reproducibility
set.seed(12342)

# Prepare spatial coordinates and basis functions
coord.grid <- down_coord_grid
coord.grid <- apply(coord.grid / 60000, 2, function(x) { x - median(x) })

basis.centers <- as.matrix(expand.grid(quantile(coord.grid[, 1], c(0.2, 0.8)), quantile(coord.grid[, 2], c(0.2, 0.8))))
basis <- lapply(1:nrow(basis.centers), function(i) {
  y <- dnorm(sqrt((coord.grid[, 1] - basis.centers[i, 1])^2 + (coord.grid[, 2] - basis.centers[i, 2])^2), mean = 0, sd = ncol(coord.grid) * 2)
  y <- y - mean(y)
  y / sqrt(sum(y^2))
})
basis <- matrix(unlist(basis), nrow = nrow(coord.grid), byrow = FALSE)
ncores <- detectCores()
unit = 2/0.03128403

# Define thresholds and modes for the analysis
thresholds_to_test <- c(0.9995)
modes_to_test <- 0:2

# --------------------------------------------------------------------------
## Jackknife Function Definition
# --------------------------------------------------------------------------

#' Perform Sequential Jackknife Standard Error Estimation for fit_model_gaus
#'
#' This function automates a sequential leave-one-out jackknife procedure.
#' Each iteration is run one by one, using the specified number of cores for the model fit.
#'
#' @param data The full input dataset vector.
#' @param coords The coordinates grid.
#' @param mode The fitting mode (0, 1, or 2).
#' @param general_init A numeric vector with all three initial parameters.
#' @param general_lb A numeric vector with all three lower bounds.
#' @param general_ub A numeric vector with all three upper bounds.
#' @param vario_func The variogram function to be used.
#' @param basis The basis functions.
#' @param full_fit_result Optional. A pre-computed fit result on the full dataset.
#' @param ncores Number of CPU cores to use for each individual model fit.
#' @param optimiser The optimisation algorithm to use.
#' @param save_path Directory to save intermediate jackknife results.
#'
#' @return A list containing the full fit, jackknife SE, and the matrix of jackknife estimates.
#'
jackknife_se_gaus_sequential <- function(
    data,
    coords,
    mode,
    general_init,
    general_lb,
    general_ub,
    vario_func,
    basis,
    full_fit_result = NULL,
    ncores = ncores,
    optimiser = "L-BFGS-B",
    save_path = "results_jackknife/"
) {
  n <- length(data)
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  # Step 1: Configure parameters based on mode
  cat(sprintf("\n--- Configuring Jackknife for Mode: %d ---\n", mode))
  if (mode == 0) {
    init_params <- general_init
    lower_bounds <- general_lb
    upper_bounds <- general_ub
  } else if (mode == 1) {
    init_params <- general_init[2]
    lower_bounds <- general_lb[2]
    upper_bounds <- general_ub[2]
  } else if (mode == 2) {
    init_params <- general_init[3]
    lower_bounds <- general_lb[3]
    upper_bounds <- general_ub[3]
  } else {
    stop("Invalid mode specified. Must be 0, 1, or 2.")
  }
  
  # Step 2: Full model fit (if not provided)
  if (is.null(full_fit_result)) {
    cat("Performing full dataset fit...\n")
    full_fit <- fit_model_gaus(
      data = data, coords = coords, init = init_params, vario_func = vario_func,
      basis = basis, ncores = ncores, lb = lower_bounds, ub = upper_bounds,
      optimiser = optimiser, mode = mode, hessian = FALSE
    )
  } else {
    cat("Using pre-computed full dataset fit.\n")
    full_fit <- full_fit_result
  }
  full_theta <- as.numeric(full_fit$par)
  cat("Full dataset parameters (theta_hat):", full_theta, "\n")
  
  # Step 3: Sequential Jackknife with Progress Bar
  cat(sprintf("Starting sequential jackknife for %d iterations (using %d cores per fit)...\n", n, ncores))
  
  results <- pblapply(1:n, function(i) {
    data_loo <- data[-i]
    fit_i <- tryCatch({
      fit_model_gaus(
        data = data_loo, coords = coords, init = full_theta, vario_func = vario_func,
        basis = basis, ncores = ncores, lb = lower_bounds, ub = upper_bounds,
        optimiser = optimiser, mode = mode, hessian = FALSE
      )
    }, error = function(e) {
      message(sprintf("Jackknife iteration %d failed: %s", i, e$message))
      return(list(par = rep(NA, length(full_theta))))
    })
    jack_theta_i <- as.numeric(fit_i$par)
    file_name <- file.path(save_path, sprintf("jackknife_mode%d_iter%04d.RData", mode, i))
    save(jack_theta_i, file = file_name)
    return(jack_theta_i)
  })
  
  # Step 4: Aggregate results and compute Standard Error
  cat("\nJackknife iterations complete. Calculating standard error...\n")
  jack_thetas <- do.call(rbind, results)
  valid_rows <- complete.cases(jack_thetas)
  n_valid <- sum(valid_rows)
  if (n_valid < n) {
    cat(sprintf("Warning: %d of %d jackknife iterations failed and were excluded.\n", n - n_valid, n))
  }
  jack_thetas_valid <- jack_thetas[valid_rows, , drop = FALSE]
  jack_var <- ((n_valid - 1) / n_valid) * colSums((sweep(jack_thetas_valid, 2, colMeans(jack_thetas_valid), "-"))^2)
  jack_se <- sqrt(jack_var)
  cat("Jackknife SE:", jack_se, "\n")
  
  return(list(
    full_fit = full_fit,
    jackknife_se = jack_se,
    jack_matrix = jack_thetas
  ))
}


# --------------------------------------------------------------------------
## Main Execution Loop
# --------------------------------------------------------------------------

for (threshold in thresholds_to_test) {
  
  cat(sprintf("\n\n##################################################################\n"))
  cat(sprintf("## Processing Threshold: %.4f\n", threshold))
  cat(sprintf("##################################################################\n"))
  
  cat("Step 1: Creating dataset...\n")
  new_dataset <- create_fitting_data(data_down_raw, threshold = threshold)
  cat("Dataset created successfully.\n")
  
  # Define general parameters
  general_init <- c(0.2, 0.2, 0.2)
  general_lb <- c(0, 0.01, 0.01)
  general_ub <- c(Inf, Inf, Inf)
  
  for (mode_to_test in modes_to_test) {
    
    cat(sprintf("\n========== Running Test for Threshold: %.4f, Mode: %d ==========\n", threshold, mode_to_test))
    
    # Step 2: Configure parameters for the specific mode
    cat("Step 2: Configuring parameters...\n")
    if (mode_to_test == 0) {
      init_params <- general_init
      lower_bounds <- general_lb
      upper_bounds <- general_ub
    } else if (mode_to_test == 1) {
      init_params <- general_init[2]
      lower_bounds <- general_lb[2]
      upper_bounds <- general_ub[2]
    } else if (mode_to_test == 2) {
      init_params <- general_init[3]
      lower_bounds <- general_lb[3]
      upper_bounds <- general_ub[3]
    }
    cat("Initial parameters:", init_params, "\n")
    cat("Lower bounds:", lower_bounds, "\n")
    cat("Upper bounds:", upper_bounds, "\n\n")
    
    # Step 3: Perform the main model fit
    cat("Step 3: Starting model fitting on full data...\n")
    fit_result <- fit_model_gaus(
      data = new_dataset,
      coords = coord.grid,
      init = init_params,
      vario_func = cal_vario_one_param,
      basis = basis,
      ncores = ncores,
      lb = lower_bounds,
      ub = upper_bounds,
      optimiser = "L-BFGS-B",
      mode = mode_to_test,
      hessian = FALSE
    )
    cat("Model fitting complete.\n")
    
    # Save and print the main fit result
    filename_fit <- sprintf("results_new1/fit_result_th%.4f_mode%d.RData", threshold, mode_to_test)
    save(fit_result, file = filename_fit)
    cat("Saved main fit result to:", filename_fit, "\n\n")
    cat("--- Fit Summary ---\n")
    print(fit_result)
    
    # Step 4: Perform Sequential Jackknife SE estimation for the current fit
    cat("\nStep 4: Starting Sequential Jackknife standard error estimation...\n")
    jack_result <- jackknife_se_gaus_sequential(
      data = new_dataset,
      coords = coord.grid,
      mode = mode_to_test,
      general_init = general_init,
      general_lb = general_lb,
      general_ub = general_ub,
      vario_func = cal_vario_one_param,
      basis = basis,
      full_fit_result = fit_result, # Reuse the initial fit
      ncores = ncores,              # Pass the total number of cores
      optimiser = "L-BFGS-B",
      save_path = sprintf("results_jackknife/th%.4f_mode%d/", threshold, mode_to_test)
    )
    
    # Save and print the jackknife result
    filename_jack <- sprintf("results_new1/jackknife_result_th%.4f_mode%d.RData", threshold, mode_to_test)
    save(jack_result, file = filename_jack)
    cat("\nSaved Jackknife result to:", filename_jack, "\n\n")
    
    cat("--- Jackknife Summary ---\n")
    cat("Full Parameters:", jack_result$full_fit$par, "\n")
    cat("Jackknife SE:   ", jack_result$jackknife_se, "\n")
    cat("==========================================================\n")
  }
}

cat("\n\nAll tasks completed.\n")