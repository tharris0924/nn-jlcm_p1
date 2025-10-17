# After your parallel processing loop completes, add this analysis code for k=3:

# Function to analyze results for both methods (k=3 version)
analyze_results_k3 <- function(method_results, method_name, n, c_rate) {
  
  # Filter out NULL results and extract successful iterations
  valid_results <- method_results[!sapply(method_results, is.null)]
  
  if(length(valid_results) == 0) {
    cat("No valid results for", method_name, "\n")
    return(NULL)
  }
  
  # Extract parameters from successful iterations
  B <- length(valid_results)
  
  # Initialize storage (k=3 dimensions)
  sigma_est <- numeric(B)
  fixef_est <- matrix(0, ncol=6, nrow=B)  # 3 components × 2 parameters each
  scale_est <- matrix(0, ncol=3, nrow=B)  # 3 components
  shape_est <- matrix(0, ncol=3, nrow=B)  # 3 components
  D_est <- matrix(0, ncol=3, nrow=B)      # vech of 2x2 matrix
  chol_est <- matrix(0, ncol=3, nrow=B)   # lower triangular elements
  omega_est <- matrix(0, ncol=2, nrow=B)  # 2 omega parameters
  CE <- numeric(B)
  MSE_y_in <- numeric(B)
  MSE_y_out <- numeric(B)
  ISE_in <- numeric(B)
  ISE_out <- numeric(B)
  IBS_in <- numeric(B)
  IBS_out <- numeric(B)
  pij <- vector("list", B)
  pij_true <- vector("list", B)
  
  # Fill arrays from results
  for(i in 1:B) {
    result <- valid_results[[i]]
    sigma_est[i] <- result$sigma_est
    fixef_est[i,] <- result$fixef_est
    scale_est[i,] <- result$scale_est
    shape_est[i,] <- result$shape_est
    if(is.matrix(result$D_est)) {
      D_est[i,] <- result$D_est
    } else {
      D_est[i,] <- result$D_est[1:3]  # Take first 3 elements
    }
    if(!is.null(result$chol_est)) {
      chol_est[i,] <- result$chol_est
    }
    if(!is.null(result$omega_est)) {
      omega_est[i,] <- result$omega_est
    }
    CE[i] <- result$CE
    MSE_y_in[i] <- result$MSE_y_in
    MSE_y_out[i] <- result$MSE_y_out
    ISE_in[i] <- result$ISE_in
    ISE_out[i] <- result$ISE_out
    IBS_in[i] <- result$IBS_in
    IBS_out[i] <- result$IBS_out
    pij[[i]] <- result$pij
    if(!is.null(result$pij_true)) {
      pij_true[[i]] <- result$pij_true
    }
  }
  
  # Get true parameters (load theta_orig_k3.RData)
  load("theta_orig_k3.RData")
  theta_orig$betas <- as.vector(t(theta_orig$betas))  # Flatten to vector
  
  # Filter out problematic estimates
  if(method_name == "NNEM") {
    idxs_to_keep <- which(sigma_est > 0)
  } else {
    idxs_to_keep <- which(sigma_est > 0)
  }
  
  # Apply filtering
  sigma_est <- sigma_est[idxs_to_keep]
  D_est <- D_est[idxs_to_keep,]
  chol_est <- chol_est[idxs_to_keep,]
  omega_est <- omega_est[idxs_to_keep,]
  fixef_est <- fixef_est[idxs_to_keep,]
  scale_est <- scale_est[idxs_to_keep,]
  shape_est <- shape_est[idxs_to_keep,]
  CE <- CE[idxs_to_keep]
  MSE_y_in <- MSE_y_in[idxs_to_keep]
  MSE_y_out <- MSE_y_out[idxs_to_keep]
  ISE_in <- ISE_in[idxs_to_keep]
  ISE_out <- ISE_out[idxs_to_keep]
  IBS_in <- IBS_in[idxs_to_keep]
  IBS_out <- IBS_out[idxs_to_keep]
  pij <- pij[idxs_to_keep]
  pij_true <- pij_true[idxs_to_keep]
  
  B <- length(idxs_to_keep)  # Update B after filtering
  
  # Calculate mixing proportions bias and MSE (for 3 components)
  true_p <- MSE_p <- numeric(3)
  
  for(k in seq_len(3)){
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
  
  # ISE calculations
  idxs_to_keep_ISE <- which(ISE_out > 0 & !is.na(ISE_out))
  ISE_out_mean <- mean(ISE_out[idxs_to_keep_ISE], na.rm = TRUE)
  idxs_to_keep_ISE_in <- which(ISE_in > 0 & !is.na(ISE_in))
  ISE_in_mean <- mean(ISE_in[idxs_to_keep_ISE_in], na.rm = TRUE)
  
  # IBS calculations
  idxs_to_keep_IBS <- which(IBS_out > 0 & !is.na(IBS_out))
  IBS_out_mean <- mean(IBS_out[idxs_to_keep_IBS], na.rm = TRUE)
  idxs_to_keep_IBS_in <- which(IBS_in > 0 & !is.na(IBS_in))
  IBS_in_mean <- mean(IBS_in[idxs_to_keep_IBS_in], na.rm = TRUE)
  
  # Create GOF matrix
  GOF <- rbind(class_err, ISE_in_mean, ISE_out_mean, IBS_in_mean, IBS_out_mean,
               mean(MSE_y_in, na.rm = TRUE), mean(MSE_y_out, na.rm = TRUE))
  rownames(GOF) <- c("CE", "ISE_in", "ISE_out","IBS_in", "IBS_out", "MSE_y_in", "MSE_y_out")
  
  # Calculate bias and MSE for parameters
  if(method_name == "NNEM") {
    # For NNEM, some parameters need square root transformation
    bias_sigma <- abs(mean(sigma_est) - theta_orig$sigma_e)
    matE <- rep(theta_orig$sigma_e, B)
    MSE_E <- Metrics::mse(matE, sigma_est)
  } else {
    # For JLCMM
    bias_sigma <- abs(mean(sigma_est) - theta_orig$sigma_e)
    matE <- rep(theta_orig$sigma_e, B)
    MSE_E <- Metrics::mse(matE, sigma_est)
  }
  
  # Fixed effects bias and MSE
  bias_fixef <- abs(colMeans(fixef_est) - matrix(theta_orig$betas, nrow = 1, byrow = T))
  matB <- matrix(rep(theta_orig$betas, B), nrow=B, ncol=6, byrow = T)
  MSE_beta <- sapply(1:6, function(i) Metrics::mse(matB[,i], fixef_est[,i]))
  
  # Variance-covariance matrix bias and MSE
  if(method_name == "NNEM") {
    bias_D <- abs(colMeans(D_est) - matrix(matrixcalc::vech(theta_orig$D), nrow = 1, byrow = T))
    matD <- matrix(rep(matrixcalc::vech(theta_orig$D), B), nrow = B, ncol=3, byrow = T)
    MSE_D <- sapply(1:3, function(i) Metrics::mse(matD[,i], D_est[,i]))
  } else {
    # For JLCMM, need to handle omega scaling
    ave_omegas <- colMeans(sqrt(omega_est^2))
    bias_D <- abs(ave_omegas[1]^2 * colMeans(D_est) - matrix(matrixcalc::vech(theta_orig$D), nrow = 1, byrow = T))
    matD <- matrix(rep(matrixcalc::vech(theta_orig$D), B), nrow = B, ncol=3, byrow = T)
    MSE_D <- sapply(1:3, function(i) Metrics::mse(matD[,i], ave_omegas[1]^2 * D_est[,i]))
  }
  
    bias_chol <- abs(colMeans(chol_est) - matrix(matrixcalc::vech(theta_orig$cholesky), nrow = 1, byrow = T))
    matchol <- matrix(rep(matrixcalc::vech(theta_orig$cholesky), B), nrow = B, ncol=3, byrow = T)
    MSE_chol <- sapply(1:3, function(i) Metrics::mse(matchol[,i], chol_est[,i]))

  
  # Weibull parameters bias and MSE
  bias_scale <- abs(colMeans(scale_est) - matrix(theta_orig$scale, nrow = 1, byrow = T))
  bias_shape <- abs(colMeans(shape_est) - matrix(theta_orig$shape, nrow = 1, byrow = T))
  
  matK <- matrix(rep(theta_orig$shape, B), nrow=B, ncol=3, byrow = T)
  MSE_K <- sapply(1:3, function(i) Metrics::mse(matK[,i], shape_est[,i]))
  
  matL <- matrix(rep(theta_orig$scale, B), nrow=B, ncol=3, byrow = T)
  MSE_L <- sapply(1:3, function(i) Metrics::mse(matL[,i], scale_est[,i]))
  
  # Omega parameters (for methods that have them)
  bias_omegas <- MSE_O <- numeric(2)
  if(method_name == "NNEM") {
    bias_omegas <- abs(colMeans(omega_est) - matrix(theta_orig$omegas, nrow = 1, byrow = T))
    matO <- matrix(rep(theta_orig$omegas, B), nrow=B, ncol=2, byrow = T)
    MSE_O <- sapply(1:2, function(i) Metrics::mse(matO[,i], omega_est[,i]))
  } else {
    # For JLCMM, omega is transformed differently
    ave_omegas <- 1/colMeans(sqrt(omega_est^2))
    bias_omegas <- abs(ave_omegas - matrix(theta_orig$omegas, nrow = 1, byrow = T))
    matO <- matrix(rep(theta_orig$omegas, B), nrow=B, ncol=2, byrow = T)
    MSE_O <- sapply(1:2, function(i) Metrics::mse(matO[,i], 1/sqrt(omega_est[,i]^2)))
  }
  
  # Create results table (k=3 version)
  results_table <- data.frame(
    Parameter = c(
      "$\\pi_1(x)$", "$\\pi_2(x)$", "$\\pi_3(x)$",
      "$\\sigma_e$",
      "$\\beta_{11}$", "$\\beta_{12}$", "$\\beta_{21}$", "$\\beta_{22}$", "$\\beta_{31}$", "$\\beta_{32}$",
      "$d_{11}$", "$d_{21}$", "$d_{22}$",
      "$\\kappa_1$", "$\\kappa_2$", "$\\kappa_3$",
      "$\\lambda_1$", "$\\lambda_2$", "$\\lambda_3$",
      "$\\omega_1$", "$\\omega_2$"
    ),
    
    True_Value = c(
      numeric(3), # Mixing proportions placeholder
      theta_orig$sigma_e,
      theta_orig$betas,
      matrixcalc::vech(theta_orig$cholesky),
      theta_orig$shape,
      theta_orig$scale,
      theta_orig$omegas
    ),
    
    Estimate = c(
      numeric(3), # Mixing proportions placeholder
      mean(sigma_est),
      colMeans(fixef_est),
      colMeans(chol_est),
      colMeans(shape_est),
      colMeans(scale_est),
      if(method_name == "NNEM") colMeans(omega_est) else 1/colMeans(sqrt(omega_est^2))
    ),
    
    Bias = c(
      # Mixing proportions bias
      mix.prop.bias,
      bias_sigma,
      bias_fixef,
      bias_chol,
      bias_shape,
      bias_scale,
      bias_omegas
    ),
    
    MSE = c(
      # Mixing proportions MSE
      MSE_p,
      MSE_E,
      MSE_beta,
      MSE_chol,
      MSE_K,
      MSE_L,
      MSE_O
    ),
    
    RMSE = c(
      # Root MSE
      sqrt(MSE_p),
      sqrt(MSE_E),
      sqrt(MSE_beta),
      sqrt(MSE_chol),
      sqrt(MSE_K),
      sqrt(MSE_L),
      sqrt(MSE_O)
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
      IBS_in = IBS_in_mean,
      IBS_out = IBS_out_mean,
      mse_y_in = mean(MSE_y_in, na.rm = TRUE),
      mse_y_out = mean(MSE_y_out, na.rm = TRUE)
    )
  ))
}

# After your foreach loop completes, analyze results for each parameter combination
library(Metrics)
library(xtable)
library(matrixcalc)

# Storage for combined results across all parameter combinations
all_nnem_results <- list()
all_jlcmm_results <- list()
combined_performance_table <- data.frame()

sample_sizes <- c(100, 500, 1000)
censoring_rates <- c(0.05, 0.25, 0.5)
combination_idx <- 0

for(n in sample_sizes) {
  for(c_rate in censoring_rates[3]) {
    combination_idx <- combination_idx + 1
    
    # Convert censoring rate to string format for file paths
    c_str <- if(c_rate == 0.05) "005" else if(c_rate == 0.25) "025" else "05"
    # NNEM\Example_k-3\jlcmm_k3_results_n100_c05.RData
    # Load the saved results for this combination
    nnem_file <- paste0("NNEM/Example_k-3/nnem_k3_results_n", n, "_c", c_str, ".RData")
    jlcmm_file <- paste0("NNEM/Example_k-3/jlcmm_k3_results_n", n, "_c", c_str, ".RData")
    
    if(file.exists(nnem_file)) {
      load(nnem_file)
    } else {
      cat("NNEM file not found:", nnem_file, "\n")
      next
    }
    
    if(file.exists(jlcmm_file)) {
      load(jlcmm_file)
    } else {
      cat("JLCMM file not found:", jlcmm_file, "\n")
      next
    }
    
    # Extract NNEM results (convert back to list format)
    nnem_results <- list()
    for(i in 1:1000) {
      if(!is.na(sigma_est[i]) && sigma_est[i] > 0) {
        nnem_results[[length(nnem_results) + 1]] <- list(
          iter = i,
          sigma_est = sigma_est[i],
          fixef_est = fixef_est[i,c(1,3,2,4,5,6)],
          scale_est = scale_est[i,],
          shape_est = shape_est[i,],
          D_est = D_est[i,],
          chol_est = chol_est[i,],
          omega_est = omega_est[i,],
          CE = CE[i],
          MSE_y_in = MSE_y_in[i],
          MSE_y_out = MSE_y_out[i],
          IBS_in=ibs_est_lcmm_in[i],
          IBS_out=ibs_est_lcmm_out[i],
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
          D_est = D_est_lcmm[i,],
          chol_est =t(chol((omega_est_lcmm[i,1]^2)*vec2sym(chol_est_lcmm[i,])))[lower.tri(matrix(0,2,2), diag=T)],
          omega_est = omega_est_lcmm[i,],
          CE = CE_lcmm[i],
          MSE_y_in = MSE_y_in_lcmm[i],
          MSE_y_out = MSE_y_out_lcmm[i],
          IBS_in=ibs_est_lcmm_in[i],
          IBS_out=ibs_est_lcmm_out[i],
          ISE_in = ise_lcmm_in[i],
          ISE_out = ise_lcmm_out[i],
          pij = pij_lcmm[[i]],
          pij_true = if(exists("pij_true") && length(pij_true) >= i) pij_true[[i]] else NULL
        )
      }
    }
    
    # Analyze results
    nnem_analysis <- analyze_results_k3(method_results = nnem_results, method_name = "NNEM", n, c_rate)
    jlcmm_analysis <- analyze_results_k3(method_results = jlcmm_results, method_name = "JLCMM", n, c_rate)
    
    # Store results
    all_nnem_results[[paste0("n", n, "_c", c_str)]] <- nnem_analysis
    all_jlcmm_results[[paste0("n", n, "_c", c_str)]] <- jlcmm_analysis
    
    # Create performance comparison table row
    if(!is.null(nnem_analysis) && !is.null(jlcmm_analysis)) {
      performance_row <- data.frame(
        Sample_Size = n,
        Censoring = c_rate,
        Method = c("NNEM", "JLCMM"),
        CE = c(nnem_analysis$performance_metrics$class_err,
               jlcmm_analysis$performance_metrics$class_err),
        ISE_In = c(nnem_analysis$performance_metrics$ISE_in,
                   jlcmm_analysis$performance_metrics$ISE_in),
        ISE_Out = c(nnem_analysis$performance_metrics$ISE_out,
                    jlcmm_analysis$performance_metrics$ISE_out),
        MSE_Y_In = c(nnem_analysis$performance_metrics$mse_y_in,
                     jlcmm_analysis$performance_metrics$mse_y_in),
        MSE_Y_Out = c(nnem_analysis$performance_metrics$mse_y_out,
                      jlcmm_analysis$performance_metrics$mse_y_out)
      )
      
      combined_performance_table <- rbind(combined_performance_table, performance_row)
    }
    
    # Generate LaTeX tables
    if(!is.null(nnem_analysis)) {
      xtable_nnem <- xtable(nnem_analysis$results_table,
                           caption = paste0("NNEM Parameter Estimates k=3 (n=", n, ", c=", c_rate, ")"),
                           label = paste0("tab:nnem_k3_n", n, "_c", c_str),
                           digits = 3)
      
      cat("\n=== NNEM TABLE k=3 n=", n, ", c=", c_rate, " ===\n")
      print(xtable_nnem, type = "latex", include.rownames = FALSE,
            sanitize.text.function = function(x) x, caption.placement = "top")
    }
    
    if(!is.null(jlcmm_analysis)) {
      xtable_jlcmm <- xtable(jlcmm_analysis$results_table,
                             caption = paste0("JLCMM Parameter Estimates k=3 (n=", n, ", c=", c_rate, ")"),
                             label = paste0("tab:jlcmm_k3_n", n, "_c", c_str),
                             digits = 3)
      
      cat("\n=== JLCMM TABLE k=3 n=", n, ", c=", c_rate, " ===\n")
      print(xtable_jlcmm, type = "latex", include.rownames = FALSE,
            sanitize.text.function = function(x) x, caption.placement = "top")
    }
  }
}

# Create combined performance comparison table
if(nrow(combined_performance_table) > 0) {
  cat("\n=== COMBINED PERFORMANCE COMPARISON k=3 ===\n")
  performance_xtable <- xtable(combined_performance_table,
                             caption = "Performance Metrics Comparison k=3: NNEM vs JLCMM",
                             label = "tab:performance_comparison_k3",
                             digits = 4)
  
  print(performance_xtable, type = "latex", include.rownames = FALSE,
        caption.placement = "top", table.placement = "H")
}

cat("Analysis completed for k=3!\n")