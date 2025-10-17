# After your parallel processing loop completes, add this analysis code:

# Function to analyze results for both methods
analyze_results <- function(method_results, method_name, n, c_rate) {
  
  # Filter out NULL results and extract successful iterations
  valid_results <- method_results[!sapply(method_results, is.null)]
  
  if(length(valid_results) == 0) {
    cat("No valid results for", method_name, "\n")
    return(NULL)
  }
  
  # Extract parameters from successful iterations
  B <- length(valid_results)
  
  # Initialize storage
  sigma_est <- numeric(B)
  fixef_est <- matrix(0, ncol=4, nrow=B)
  scale_est <- matrix(0, ncol=2, nrow=B)
  shape_est <- matrix(0, ncol=2, nrow=B)
  D_est <- numeric(B)
  chol_est <- numeric(B)
  CE <- numeric(B)
  MSE_y_in <- numeric(B)
  MSE_y_out <- numeric(B)
  ISE_in <- numeric(B)
  ISE_out <- numeric(B)
  pij <- vector("list", B)
  pij_true <- vector("list", B)
  
  # Fill arrays from results
  for(i in 1:B) {
    result <- valid_results[[i]]
    sigma_est[i] <- result$sigma_est
    fixef_est[i,] <- result$fixef_est
    scale_est[i,] <- result$scale_est
    shape_est[i,] <- result$shape_est
    D_est[i] <- result$D_est
    chol_est[i] <- result$chol_est
    CE[i] <- result$CE
    MSE_y_in[i] <- result$MSE_y_in
    MSE_y_out[i] <- result$MSE_y_out
    ISE_in[i] <- result$ISE_in
    ISE_out[i] <- result$ISE_out
    pij[[i]] <- result$pij
    if(!is.null(result$pij_true)) {
      pij_true[[i]] <- result$pij_true
    }
  }
  
  # Get true parameters (assuming they're the same across iterations)
  # You'll need to load one result file to get theta_orig
  load(paste0("NNEM/Example2/k=2_Results_n=", n, "_c", c_rate, "/nnem_size1", n, ".RData"))
  theta_orig <- result$theta_orig
  theta_orig$betas <- c(t(theta_orig$betas))
  
  # Filter out problematic estimates (following your original logic)
  if(method_name == "NNEM") {
    idxs_to_keep <- which(sigma_est > 0)
  } else {
    idxs_to_keep <- which(sigma_est > 0)
  }
  
  # Apply filtering
  sigma_est <- sigma_est[idxs_to_keep]
  D_est <- D_est[idxs_to_keep]
  fixef_est <- fixef_est[idxs_to_keep,]
  scale_est <- scale_est[idxs_to_keep,]
  shape_est <- shape_est[idxs_to_keep,]
  chol_est <- chol_est[idxs_to_keep,]
  CE <- CE[idxs_to_keep]
  MSE_y_in <- MSE_y_in[idxs_to_keep]
  MSE_y_out <- MSE_y_out[idxs_to_keep]
  ISE_in <- ISE_in[idxs_to_keep]
  ISE_out <- ISE_out[idxs_to_keep]
  
  # Calculate mixing proportions bias and MSE
  true_p <- MSE_p <- numeric(2)
  
  for(k in seq_len(2)){
    MAB_pi <- true_pi <- matrix(ncol=B, nrow = n)
    
    for(i in seq_len(B)){
      if(!is.null(pij[[i]]) && !is.null(pij_true[[i]])){
        # Check if pij has the correct dimensions
        if(is.matrix(pij[[i]]) && ncol(pij[[i]]) >= k) {
          MAB_pi[,i] <- pij[[i]][,k]
        }
        if(is.matrix(pij_true[[i]]) && ncol(pij_true[[i]]) >= k) {
          true_pi[,i] <- pij_true[[i]][,k]
        }
      }
    }
    
    # Calculate bias and MSE for mixing proportions
    pi2 <- mean(abs(rowMeans(MAB_pi, na.rm = T) - rowMeans(true_pi, na.rm = T)), na.rm = T)
    pi3 <- mean(rowMeans(((MAB_pi) - (true_pi))^2, na.rm = T), na.rm = T)
    true_p[k] <- pi2
    MSE_p[k] <- pi3
  }
  
  mix.prop.bias <- true_p
  
  # Calculate performance metrics
  class_err <- mean(CE, na.rm = TRUE)
  
  # IBS calculations
  idxs_to_keep_ISE <- which(ISE_out > 0)
  ISE_out_mean <- mean(ISE_out[idxs_to_keep_ISE], na.rm = TRUE)
  idxs_to_keep_ISE_in <- which(ISE_in > 0)
  ISE_in_mean <- mean(ISE_in[idxs_to_keep_ISE_in], na.rm = TRUE)
  
  # Create GOF matrix
  GOF <- rbind(class_err, ISE_in_mean, ISE_out_mean, 0, 0) # Last two are placeholders for ISE
  rownames(GOF) <- c("CE", "ISE_in", "ISE_out", "ISE_in", "ISE_out")
  
  # Print summary statistics
  # cat("\n", method_name, "Results for n =", n, ", c_rate =", c_rate, "\n")
  # cat("Valid iterations:", B, "out of 1000\n")
  # cat("Error Variance:", mean(sqrt(sigma_est)), "vs", theta_orig$sigma_e, "\n")
  # cat("RE Variance-Covariance:", mean(sqrt(D_est)), "vs", sqrt(theta_orig$D[1]), "\n")
  # cat("Longitudinal Fixef:", colMeans(fixef_est), "vs", theta_orig$betas, "\n")
  # cat("Weibull Scales:", colMeans(scale_est), "vs", theta_orig$scale, "\n")
  # cat("Weibull Shapes:", colMeans(shape_est), "vs", theta_orig$shape, "\n")
  
  # Calculate bias and MSE
  if(method_name == "NNEM") {
    bias_sigma <- abs(mean(sqrt(sigma_est)) - theta_orig$sigma_e)
    matE <- rep(theta_orig$sigma_e, B)
    MSE_E <- Metrics::mse(matE, sqrt(sigma_est))
  } else {
    bias_sigma <- abs(mean(sigma_est) - theta_orig$sigma_e)
    matE <- rep(theta_orig$sigma_e, B)
    MSE_E <- Metrics::mse(matE, sigma_est)
  }
  
  bias_scale <- abs(colMeans(scale_est) - matrix(theta_orig$scale, nrow = 1, byrow = T))
  bias_shape <- abs(colMeans(shape_est) - matrix(theta_orig$shape, nrow = 1, byrow = T))
  bias_fixef <- abs(colMeans(fixef_est) - matrix(theta_orig$betas, nrow = 1, byrow = T))
  
  if(method_name == "NNEM") {
    bias_D <- abs(mean(sqrt(D_est)) - sqrt(theta_orig$D[1]))
    matD <- matrix(theta_orig$D, nrow = B, ncol=1, byrow = T)
    MSE_D <- Metrics::mse(sqrt(matD), sqrt(D_est))
  } else {
    bias_D <- abs(mean(D_est) - sqrt(theta_orig$D[1]))
    matD <- matrix(theta_orig$D, nrow = B, ncol=1, byrow = T)
    MSE_D <- Metrics::mse(sqrt(matD), D_est)
  }

    if(method_name == "NNEM") {
    bias_D <- abs(mean(sqrt(D_est)) - sqrt(theta_orig$D[1]))
    mat_chol <- matrix(theta_orig$cholesky[lower.tri(theta_orig$cholesky, diag=T)], nrow = B, ncol=1, byrow = T)
    MSE_D <- Metrics::mse(sqrt(matD), sqrt(D_est))
  } else {
    bias_D <- abs(mean(D_est) - sqrt(theta_orig$D[1]))
    matD <- matrix(theta_orig$D, nrow = B, ncol=1, byrow = T)
    MSE_D <- Metrics::mse(sqrt(matD), D_est)
  }
  
  # Calculate MSE for other parameters
  matB <- matrix(rep(theta_orig$betas, B), nrow=B, ncol=4, byrow = T)
  MSE_beta <- sapply(1:4, function(i) Metrics::mse(matB[,i], fixef_est[,i]))
  
  matK <- matrix(rep(theta_orig$shape, B), nrow=B, ncol=2, byrow = T)
  MSE_K <- sapply(1:2, function(i) Metrics::mse(matK[,i], shape_est[,i]))
  
  matL <- matrix(rep(theta_orig$scale, B), nrow=B, ncol=2, byrow = T)
  MSE_L <- sapply(1:2, function(i) Metrics::mse(matL[,i], scale_est[,i]))
  
  # Create results table
  results_table <- data.frame(
    Parameter = c(
      "$\\pi_1(x)$", "$\\pi_2(x)$",
      "$\\sigma_e$",
      "$\\beta_{11}$", "$\\beta_{12}$", "$\\beta_{21}$", "$\\beta_{22}$",
      "$d_{11}$",
      "$\\kappa_1$", "$\\kappa_2$",
      "$\\lambda_1$", "$\\lambda_2$"
    ),
    
    True_Value = c(
      numeric(2), # Mixing proportions placeholder
      theta_orig$sigma_e,
      theta_orig$betas,
      sqrt(theta_orig$D[1]),
      theta_orig$shape,
      theta_orig$scale
    ),
    
    Estimate = c(
      numeric(2), # Mixing proportions placeholder
      if(method_name == "NNEM") mean(sqrt(sigma_est)) else mean(sigma_est),
      colMeans(fixef_est),
      if(method_name == "NNEM") mean(sqrt(D_est)) else mean(D_est),
      colMeans(shape_est),
      colMeans(scale_est)
    ),
    
    Bias = c(
      # Mixing proportions bias
      mix.prop.bias,
      bias_sigma,
      bias_fixef,
      bias_D,
      bias_shape,
      bias_scale
    ),
    
    MSE = c(
      # Mixing proportions MSE
      MSE_p,
      MSE_E,
      MSE_beta,
      MSE_D,
      MSE_K,
      MSE_L
    ),
    
    RMSE = c(
      # Root MSE for mixing proportions
      sqrt(MSE_p),
      sqrt(MSE_E),
      sqrt(MSE_beta),
      sqrt(MSE_D),
      sqrt(MSE_K),
      sqrt(MSE_L)
    )
  )
  
  # Round numeric columns
  numeric_cols <- c("True_Value", "Estimate", "Bias", "MSE", "RMSE")
  results_table[numeric_cols] <- lapply(results_table[numeric_cols], function(x) round(x, 3))
  results_table$True_Value <- sprintf("%.10g", results_table$True_Value)
  
  return(list(
    results_table = results_table,
    GOF = GOF,
    performance_metrics = list(
      class_err = class_err,
      ISE_in = ISE_in_mean,
      ISE_out = ISE_out_mean,
      mse_y_in = mean(MSE_y_in, na.rm = TRUE),
      mse_y_out = mean(MSE_y_out, na.rm = TRUE)
    )
  ))
}

# After your foreach loop completes, analyze results for each parameter combination
library(Metrics)
library(xtable)

# Storage for combined results across all parameter combinations
all_nnem_results <- list()
all_jlcmm_results <- list()
combined_performance_table <- data.frame()
sample_sizes <- c(100,500,1000)
censoring_rates <- c(0.05,0.25,0.5)
combination_idx <- 0
for(c_rate in censoring_rates) {
for(n in sample_sizes) {

    combination_idx <- combination_idx + 1
    
    # Load the saved results for this combination
    c_str <- if(c_rate == 0.05) "005" else if(c_rate == 0.25) "025" else "05"

    load(paste0("NNEM/Example2/nnem_results_n", n, "_c", c_str, ".RData"))
    load(paste0("NNEM/Example2/jlcmm_results_n", n, "_c", c_str, ".RData"))
    
    # Extract NNEM results (convert back to list format)
    nnem_results <- list()
    for(i in 1:1000) {
      if(!is.na(sigma_est[i]) && sigma_est[i] > 0) {
        nnem_results[[length(nnem_results) + 1]] <- list(
          iter = i,
          sigma_est = sigma_est[i],
          fixef_est = fixef_est[i,],
          scale_est = scale_est[i,],
          shape_est = shape_est[i,],
          D_est = D_est[i],
          CE = CE[i],
          MSE_y_in = MSE_y_in[i],
          MSE_y_out = MSE_y_out[i],
          ISE_in = ISE_in[i],
          ISE_out = ISE_out[i],
          pij = pij[[i]],
          pij_true = pij_true[[i]]
        )
      }
    }
    
    # Extract JLCMM results (convert back to list format)
    jlcmm_results <- list()
    for(i in 1:1000) {
      if(!is.na(sigma_est_lcmm[i]) && sigma_est_lcmm[i] > 0) {
        jlcmm_results[[length(jlcmm_results) + 1]] <- list(
          iter = i,
          sigma_est = sigma_est_lcmm[i],
          fixef_est = fixef_est_lcmm[i,],
          scale_est = scale_est_lcmm[i,],
          shape_est = shape_est_lcmm[i,],
          D_est = sqrt(D_est_lcmm[i]),
          CE = CE_lcmm[i],
          MSE_y_in = MSE_y_in_lcmm[i],
          MSE_y_out = MSE_y_out_lcmm[i],
          ISE_in = ISE_in_lcmm[i],
          ISE_out = ISE_out_lcmm[i],
          pij = pij_lcmm[[i]],
          pij_true = if(exists("pij_true") && length(pij_true) >= i) pij_true[[i]] else NULL
        )
      }
    }
    
    # Analyze results
    nnem_analysis <- analyze_results(nnem_results, "NNEM", n, c_rate)
    jlcmm_analysis <- analyze_results(jlcmm_results, "JLCMM", n, c_rate)
    
    # Store results
    all_nnem_results[[paste0("n", n, "_c", c_str)]] <- nnem_analysis
    all_jlcmm_results[[paste0("n", n, "_c", c_str)]] <- jlcmm_analysis
    
    # Create performance comparison table row
    performance_row <- data.frame(
      Sample_Size = n,
      Censoring = c_rate,
      Method = c("NNEM", "JLCMM"),
      CE = c(nnem_analysis$performance_metrics$class_err,
             jlcmm_analysis$performance_metrics$class_err),
      IBS_In = c(nnem_analysis$performance_metrics$ISE_in,
                 jlcmm_analysis$performance_metrics$ISE_in),
      IBS_Out = c(nnem_analysis$performance_metrics$ISE_out,
                  jlcmm_analysis$performance_metrics$ISE_out),
      MSE_Y_In = c(nnem_analysis$performance_metrics$mse_y_in,
                   jlcmm_analysis$performance_metrics$mse_y_in),
      MSE_Y_Out = c(nnem_analysis$performance_metrics$mse_y_out,
                    jlcmm_analysis$performance_metrics$mse_y_out)
    )
    
    combined_performance_table <- rbind(combined_performance_table, performance_row)
    
    # Generate LaTeX tables
    if(!is.null(nnem_analysis)) {
      xtable_nnem <- xtable(nnem_analysis$results_table,
                           caption = paste0("NNEM Parameter Estimates (n=", n, ", c=", c_rate, ")"),
                           label = paste0("tab:nnem_n", n, "_c", c_str),
                           digits = 3)
      
      cat("\n=== NNEM TABLE n=", n, ", c=", c_rate, " ===\n")
      print(xtable_nnem, type = "latex", include.rownames = FALSE,
            sanitize.text.function = function(x) x, caption.placement = "top")
    }
    
    if(!is.null(jlcmm_analysis)) {
      xtable_jlcmm <- xtable(jlcmm_analysis$results_table,
                             caption = paste0("JLCMM Parameter Estimates (n=", n, ", c=", c_rate, ")"),
                             label = paste0("tab:jlcmm_n", n, "_c", c_str),
                             digits = 3)
      
      cat("\n=== JLCMM TABLE n=", n, ", c=", c_rate, " ===\n")
      print(xtable_jlcmm, type = "latex", include.rownames = FALSE,
            sanitize.text.function = function(x) x, caption.placement = "top")
    }
  }
}

# Create combined performance comparison table
cat("\n=== COMBINED PERFORMANCE COMPARISON ===\n")
performance_xtable <- xtable(combined_performance_table,
                           caption = "Performance Metrics Comparison: NNEM vs JLCMM",
                           label = "tab:performance_comparison",
                           digits = 4)

print(performance_xtable, type = "latex", include.rownames = FALSE,
      caption.placement = "top", table.placement = "H")
