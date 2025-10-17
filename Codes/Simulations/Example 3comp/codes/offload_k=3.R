# Define parameter combinations for k=3
sample_sizes <- c(100, 500, 1000)
censoring_rates <- c(0.05, 0.25, 0.5)

library(doParallel)
library(lcmm)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
source("C:/Users/th240/Documents/Ph.D. Statistics/Codes/NNEM/Example_k-3/estep.R")
source("C:/Users/th240/Documents/Ph.D. Statistics/Codes/NNEM/Example_k-3/bias_utils.R")
source("C:/Users/th240/Documents/Ph.D. Statistics/Codes/NNEM/Example_k-3/em_utils.R")
load("theta_orig_k3.RData")

# Create progress tracking
total_combinations <- length(sample_sizes) * length(censoring_rates)
current_combination <- 0
# n <- 100
# c_rate <- 0.05
# Outer loop for sample sizes
for(n in sample_sizes[1]) {
  # Inner loop for censoring rates
  for(c_rate in censoring_rates[1]) {
    current_combination <- current_combination + 1
    
    # Progress bar for parameter combinations
    cat("\n", rep("=", 60), "\n", sep="")
    cat("PARAMETER COMBINATION", current_combination, "of", total_combinations, "\n")
    cat("n =", n, ", censoring rate =", c_rate, "\n")
    progress_pct <- round((current_combination-1) / total_combinations * 100)
    progress_bar <- paste0("[", 
                           paste(rep("█", floor(progress_pct/2)), collapse=""),
                           paste(rep("░", 50 - floor(progress_pct/2)), collapse=""),
                           "] ", progress_pct, "%")
    cat(progress_bar, "\n")
    cat(rep("=", 60), "\n\n", sep="")
    
    cat("Processing n =", n, ", censoring rate =", c_rate, "\n")
    
    # Initialize result storage for this combination
    B <- 1000
    
    # NNEM results (k=3 parameters)
    sigma_est <- matrix(0, nrow = B)
    fixef_est <- matrix(0, ncol=6, nrow=B)  # 3 components × 2 parameters each
    scale_est <- matrix(0, ncol=3, nrow=B)  # 3 components
    shape_est <- matrix(0, ncol=3, nrow=B)  # 3 components
    D_est <- matrix(0, ncol=3, nrow=B)      # vech of 2x2 matrix = 3 elements
    chol_est <- matrix(0, ncol=3, nrow=B)   # lower triangular elements
    omega_est <- matrix(0, ncol=2, nrow=B)  # 2 omega parameters
    pij <- c_true <- pred_c <- list()
    pij_true <- list()
    CE <- ibs_est_train <- ibs_est_test <- ISE_in <- ISE_out <- MSE_y_out <- MSE_y_in <- numeric(B)
    
    # JLCMM results (k=3 parameters - same structure as NNEM)
    sigma_est_lcmm <- matrix(0, nrow = B)
    fixef_est_lcmm <- matrix(0, ncol=6, nrow=B)
    scale_est_lcmm <- matrix(0, ncol=3, nrow=B)
    shape_est_lcmm <- matrix(0, ncol=3, nrow=B)
    chol_est_lcmm <- matrix(0, ncol=3, nrow=B)
    D_est_lcmm <- matrix(0, ncol=3, nrow=B)
    omega_est_lcmm <- matrix(0, ncol=2, nrow=B)
    pij_lcmm <- pij_nnem <- list()
    CE_lcmm <- ibs_est_lcmm_in <- ibs_est_lcmm_out <- ise_lcmm_in <- ise_lcmm_out <- MSE_y_out_lcmm <- MSE_y_in_lcmm <- numeric(B)
    
    # Convert censoring rate to string format for file paths
    c_str <- if(c_rate == 0.05) "0.05" else if(c_rate == 0.25) "0.25" else "0.5"
    
    # Parallel processing for both NNEM and JLCMM
    offload <- foreach(
      iter = 1:1000,
      .packages = c("matrixcalc", "prodlim", "data.table", "survival", "SurvMetrics", "lcmm"),
      .export = ls(),
      .combine = c, 
      .verbose = TRUE, .errorhandling = "pass"
    ) %dopar% {
      # iter <- 1
      # Load both NNEM and JLCMM results
      wkdir_nnem <- paste0("NNEM/Example_k-3/k=3_Results_n=", n, "_c", c_str)
      file_nnem <- paste0(wkdir_nnem, "/nnem_size", iter, n, ".RData")
      file_jlcmm <- paste0(wkdir_nnem, "/jlcmm_fit", iter, n, ".RData")
      
      # Initialize result containers
      nnem_result <- NULL
      jlcmm_result <- NULL
      
      # Process NNEM results
      if(file.exists(file_nnem)) {
        load(file_nnem)
        theta <- result$theta
        
        if(!is.null(theta)){
          # Extract NNEM parameters for k=3
          local_sigma_est <- theta$sigma_e
          local_fixef_est <- matrix(theta$betas, nrow=1, byrow=T)  # 6 parameters
          local_scale_est <- matrix(theta$scale, nrow=1, byrow=T)  # 3 components
          local_shape_est <- matrix(theta$shape, nrow=1, byrow=T) # 3 components
          local_D_est <- matrixcalc::vech(theta$D)                # 3 elements
          local_chol_est <- theta$cholesky[lower.tri(theta$cholesky, diag=T)]  # 3 elements
          local_omega_est <- matrix(sqrt((theta$omegas)^2), nrow=1, byrow=T)   # 2 elements
          local_pij <- theta$pj
          
          # Get true components
          true_components <- numeric(n)
          for(j in seq_len(n)){
            red_c <- result$data_long[result$data_long$ID==j,]
            true_components[j] <- unique(red_c$g)
          }
          
          # True mixing proportions for k=3
          true_mixing_prop <- function(x) {
            pi1 <- (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4)) * as.numeric(x < 3)
            pi2 <- (1 - (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4))) * as.numeric(x < 0)
            pi3 <- (2 * exp(-0.1 * (x - 6)^4)) / (1 + exp(-0.1 * (x - 6)^4)) * as.numeric(x >= 3)
            pi4 <- 1 - (pi1 + pi3)
            return(cbind(pi1, pi3, pi4))
          }
          
          local_pij_true <- true_mixing_prop(result$x)
          
          # Classification error
          local_CE <- mean(result$C != true_components)
          
          # In-sample metrics
          data <- result$data_long
          data_wide <- result$data_wide
          comp.p <- estep(local_pij, theta)$P
          
          local_ibs_est_train <- get_nnem_brier_score(
            shape_params=local_shape_est, 
            scale_params=local_scale_est, 
            mix.prob=comp.p, 
            test_data=result$data_wide
          )$IBS$ibs
          
          true_comp.p <- estep(local_pij_true, theta_orig)$P
          local_ISE_in <- get_nnem_ISE(result, 
                                       shapes = local_shape_est, 
                                       scales = local_scale_est, 
                                       mix.prob = comp.p,
                                       theta_orig = theta_orig, 
                                       true_pi = true_comp.p)
          
          # Longitudinal prediction for NNEM
          pred_y <- matrix(0, nrow(result$data_long), 3)  # 3 components
          Z <- model.matrix(~measure_time, data = result$data_long)
          X <- model.matrix(~measure_time, data = result$data_long)
          
          n_timepoints <- table(result$data_long$ID)     
          bi <- list()
          m <- length(theta$shape)  # Should be 3
          q <- nrow(theta$D)
          
          # Compute class-specific random effects for NNEM
          for(j in 1:m){
            alphas <- c(theta$omegas, 1)
            sigma_g <- theta$sigma_e
            D_g <- (alphas[j]^2) * theta$D
            beta_g <- theta$betas[j,]
            
            bigs <- matrix(nrow = n, ncol = 2)
            for(k in 1:n){
              n_i <- n_timepoints[k]
              sim <- data
              subset <- sim[sim$ID==k,]
              Z_i <- Z[sim$ID==k, , drop = FALSE]
              X_i <- X[sim$ID==k, , drop = FALSE]
              Xi_beta <- X_i %*% beta_g
              Y_i <- subset$y
              
              Vi <- as.matrix(Z_i %*% D_g %*% t(Z_i) + sigma_g^2 * diag(n_i))
              big <- D_g %*% t(Z_i) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
              bigs[k,] <- big
            }
            bi[j] <- list(bigs)
          }
          
          # Predict longitudinal outcomes for each component
          for(g in 1:3){
            mu <- NULL
            beta_tmp <- theta$betas[g,]
            ranef_tmp <- bi[[g]]
            
            for (k in 1:n) {
              idx <- which(result$data_long$ID == k)
              Z_i <- Z[idx, , drop = FALSE]
              X_i <- X[idx, , drop = FALSE]
              XBg <- (X_i %*% beta_tmp + Z_i %*% ranef_tmp[k,])
              mu <- c(mu, XBg)
            }
            pred_y[,g] <- mu
          }
          
          preds <- rowSums(pred_y * comp.p[result$data_long$ID,])
          local_MSE_y_in <- Metrics::mse(actual = result$data_long$y, predicted = preds)
          
          # Out-of-sample evaluation for NNEM
          if(iter==1000){
            file_test <- paste0(wkdir_nnem,"/nnem_size", 1, n,".RData")
          } else {
            file_test <- paste0(wkdir_nnem,"/nnem_size", iter+1, n,".RData")
          }
          
          local_MSE_y_out <- NA
          local_ibs_est_test <- NA
          local_ISE_out <- NA
          
          if(file.exists(file_test)) {
            load(file_test)
            if(!is.null(result)){
              local_pij_true <- true_mixing_prop(result$x)
              data <- result$data_long
              data_wide <- result$data_wide
              true_comp.p <- estep(local_pij_true, theta_orig)$P
              test_comp.p <- estep(prior_mixing = theta$pj, theta)$P
              
              local_ISE_out <- get_nnem_ISE(result, 
                                            shapes = local_shape_est, 
                                            scales = local_scale_est, 
                                            mix.prob = test_comp.p,
                                            theta_orig = theta_orig, 
                                            true_pi = true_comp.p)
              
              local_ibs_est_test <- get_nnem_brier_score(
                shape_params = local_shape_est, 
                scale_params = local_scale_est, 
                mix.prob = test_comp.p, 
                test_data = result$data_wide
              )$IBS$ibs
              
              # Out-of-sample longitudinal prediction
              pred_y <- matrix(0, nrow(result$data_long), 3)
              Z <- model.matrix(~measure_time, data = result$data_long)
              X <- model.matrix(~measure_time, data = result$data_long)
              
              n_timepoints <- table(result$data_long$ID)
              bi <- list()
              
              for(j in 1:m){
                alphas <- c(theta$omegas, 1)
                sigma_g <- theta$sigma_e
                D_g <- (alphas[j]^2) * theta$D
                beta_g <- theta$betas[j,]
                
                bigs <- matrix(nrow = n, ncol = 2)
                for(k in 1:n){
                  n_i <- n_timepoints[k]
                  sim <- data
                  subset <- sim[sim$ID==k,]
                  Z_i <- Z[sim$ID==k, , drop = FALSE]
                  X_i <- X[sim$ID==k, , drop = FALSE]
                  Xi_beta <- X_i %*% beta_g
                  Y_i <- subset$y
                  
                  Vi <- as.matrix(Z_i %*% D_g %*% t(Z_i) + sigma_g^2 * diag(n_i))
                  big <- D_g %*% t(Z_i) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
                  bigs[k,] <- big
                }
                bi[j] <- list(bigs)
              }
              
              for(g in 1:3){
                mu <- NULL
                beta_tmp <- theta$betas[g,]
                ranef_tmp <- bi[[g]]
                
                for (k in 1:n) {
                  idx <- which(result$data_long$ID == k)
                  Z_i <- Z[idx, , drop = FALSE]
                  X_i <- X[idx, , drop = FALSE]
                  XBg <- (X_i %*% beta_tmp + Z_i %*% ranef_tmp[k,])
                  mu <- c(mu, XBg)
                }
                pred_y[,g] <- mu
              }
              
              preds <- rowSums(pred_y * test_comp.p[result$data_long$ID,])
              local_MSE_y_out <- Metrics::mse(actual = result$data_long$y, predicted = preds)
            }
          }
          
          nnem_result <- list(
            method = "NNEM",
            iter = iter,
            sigma_est = local_sigma_est,
            fixef_est = local_fixef_est,
            scale_est = local_scale_est,
            shape_est = local_shape_est,
            D_est = local_D_est,
            chol_est = local_chol_est,
            omega_est = local_omega_est,
            pij = local_pij,
            pij_true = local_pij_true,
            MSE_y_in = local_MSE_y_in,
            MSE_y_out = local_MSE_y_out,
            ibs_est_train = local_ibs_est_train,
            ibs_est_test = local_ibs_est_test,
            ISE_in = local_ISE_in,
            ISE_out = local_ISE_out,
            CE = local_CE
          )
        }
      }
      
      # Process JLCMM results
      if(file.exists(file_jlcmm) && file.exists(file_nnem)) {
        load(file_jlcmm)  # loads jlcmm_fit
        load(file_nnem)   # loads result
        
        if(!is.null(result) && !is.null(jlcmm_fit) && !is.na(jlcmm_fit$best[1])) {
          est <- jlcmm_fit$best
          obj <- jlcmm_fit
          nclasses <- ncol(obj$pprob) - 2
          coefs <- obj$best
          coefend <- min(which(grepl('Weibull', names(coefs)))) - 1
          coefs_multilogit <- matrix(coefs[1:coefend], nrow = nclasses - 1)
          tmpX <- cbind(1, result$x)
          linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
          exps <- exp(cbind(linearval, 0))
          pi_logit <- exps / rowSums(exps)
          
          # True components
          true_components <- numeric(n)
          for(j in seq_len(n)){
            red_c <- result$data_long[result$data_long$ID==j,]
            true_components[j] <- unique(red_c$g)
          }
          
          # True mixing proportions
          true_mixing_prop <- function(x) {
            pi1 <- (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4)) * as.numeric(x < 3)
            pi2 <- (1 - (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4))) * as.numeric(x < 0)
            pi3 <- (2 * exp(-0.1 * (x - 6)^4)) / (1 + exp(-0.1 * (x - 6)^4)) * as.numeric(x >= 3)
            pi4 <- 1 - (pi1 + pi3)
            return(cbind(pi1, pi3, pi4))
          }
          
          
          pi_true <- true_mixing_prop(result$x)
          local_CE_lcmm <- mean(jlcmm_fit$pprob$class != true_components)
          
          # Extract JLCMM parameters
          local_sigma_est_lcmm <- est["stderr"]
          local_fixef_est_lcmm <- est[c(13, 16, 12, 15, 11, 14)]  # Based on your bias script
          local_shape_est_lcmm <- exp(est[c(10, 8, 6)])
          local_scale_est_lcmm <- exp(est[c(9, 7, 5)])^(-1 / exp(est[c(10, 8, 6)]))
          
          # Variance parameters
          idxs_varcov <- which(grepl('varcov', names(est)))
          idxs_varprop <- which(grepl('varprop', names(est)))
          local_chol_est_lcmm <- est[idxs_varcov]
          local_D_est_lcmm <- coefs[idxs_varcov]
          local_omega_est_lcmm <- est[idxs_varprop]
          
          # In-sample evaluation
          # result$data_wide$x <- NULL
          local_ise_lcmm_in <- get_nnem_ISE(result, 
                                            shapes = local_shape_est_lcmm, 
                                            scales = local_scale_est_lcmm, 
                                            mix.prob = as.matrix(jlcmm_fit$pprob[,3:5]),
                                            theta_orig = theta_orig, 
                                            true_pi = true_comp.p)
          
          local_ibs_est_lcmm_in  <- get_nnem_brier_score(
            shape_params=local_shape_est_lcmm, 
            scale_params=local_scale_est_lcmm, 
            mix.prob=as.matrix(jlcmm_fit$pprob[,3:5]), 
            test_data=result$data_wide
          )$IBS$ibs
          
          # Longitudinal predictions for JLCMM
          mod <- jlcmm_fit
          predy <- mod$pred$pred_ss
          local_MSE_y_in_lcmm <- Metrics::mse(result$data_long$y, predy)
          
          # Out-of-sample evaluation for JLCMM
          if(iter==1000){
            file_test <- paste0(wkdir_nnem,"/nnem_size", 1, n,".RData")
          } else {
            file_test <- paste0(wkdir_nnem,"/nnem_size", iter+1, n,".RData")
          }
          
          local_MSE_y_out_lcmm <- NA
          local_ibs_est_lcmm_out <- NA
          local_ise_lcmm_out <- NA
          
          if(file.exists(file_test)) {
            load(file_test)
            if(!is.null(result)){
              tmpX <- cbind(1, result$x)
              linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
              exps <- exp(cbind(linearval, 0))
              pi_logit_test <- exps / rowSums(exps)
              
              pred_mixing <- estep(pi_logit_test, theta)$P
              
             ( local_ibs_est_lcmm_out <- get_nnem_brier_score(
            shape_params=local_shape_est_lcmm, 
            scale_params=local_scale_est_lcmm, 
            mix.prob= pred_mixing,
            test_data=result$data_wide
          )$IBS$ibs)
              
              local_ise_lcmm_out <- get_nnem_ISE(result, 
                                            shapes = local_shape_est_lcmm, 
                                            scales = local_scale_est_lcmm, 
                                            mix.prob = pred_mixing,
                                            theta_orig = theta_orig, 
                                            true_pi = true_comp.p)
              
              # Out-of-sample longitudinal prediction using JLCMM
     
              
               result$data_long$x <- result$x[result$data_long$ID]
              predy_test_raw <- rowSums(pi_logit_test[result$data_long$ID,] * 
                                        lcmm::predictY(jlcmm_fit, newdata=result$data_long)$pred)
              local_MSE_y_out_lcmm <- Metrics::mse(result$data_long$y, predy_test_raw)
            }
          }
          
          jlcmm_result <- list(
            method = "JLCMM",
            iter = iter,
            sigma_est = local_sigma_est_lcmm,
            fixef_est = local_fixef_est_lcmm,
            scale_est = local_scale_est_lcmm,
            shape_est = local_shape_est_lcmm,
            D_est = local_D_est_lcmm,
            chol_est = local_chol_est_lcmm,
            omega_est = local_omega_est_lcmm,
            pij = pi_logit,
            pij_nnem = result$theta$pj,
            pij_true = pi_true,
            MSE_y_in = local_MSE_y_in_lcmm,
            MSE_y_out = local_MSE_y_out_lcmm,
            ibs_est_train = local_ibs_est_lcmm_in,
            ibs_est_test = local_ibs_est_lcmm_out,
            ISE_in = local_ise_lcmm_in,
            ISE_out = local_ise_lcmm_out,
            CE = local_CE_lcmm
          )
        }
      }
      
      # Return both results
      return(list(nnem = nnem_result, jlcmm = jlcmm_result))
    }
    
    # Process results and assign to appropriate variables
    for(result_pair in offload) {
      # Process NNEM results
      if(is.list(result_pair)){
      if(!is.null(result_pair$nnem)) {
        result <- result_pair$nnem
        iter <- result$iter
        sigma_est[iter] <- result$sigma_est
        fixef_est[iter,] <- result$fixef_est
        scale_est[iter,] <- result$scale_est
        shape_est[iter,] <- result$shape_est
        D_est[iter,] <- result$D_est
        chol_est[iter,] <- result$chol_est
        omega_est[iter,] <- result$omega_est
        pij[[iter]] <- result$pij
        pij_true[[iter]] <- result$pij_true
        MSE_y_in[iter] <- result$MSE_y_in
        MSE_y_out[iter] <- result$MSE_y_out
        ibs_est_train[iter] <- result$ibs_est_train
        ibs_est_test[iter] <- result$ibs_est_test
        ISE_in[iter] <- result$ISE_in
        ISE_out[iter] <- result$ISE_out
        CE[iter] <- result$CE
      }}
      
      # Process JLCMM results
      if(is.list(result_pair)){
      if(!is.null(result_pair$jlcmm)) {
        result <- result_pair$jlcmm
        iter <- result$iter
        sigma_est_lcmm[iter] <- result$sigma_est
        fixef_est_lcmm[iter,] <- result$fixef_est
        scale_est_lcmm[iter,] <- result$scale_est
        shape_est_lcmm[iter,] <- result$shape_est
        D_est_lcmm[iter,] <- result$D_est
        chol_est_lcmm[iter,] <- result$chol_est
        omega_est_lcmm[iter,] <- result$omega_est
        pij_lcmm[[iter]] <- result$pij
        pij_nnem[[iter]] <- result$pij_nnem
        MSE_y_in_lcmm[iter] <- result$MSE_y_in
        MSE_y_out_lcmm[iter] <- result$MSE_y_out
        ibs_est_lcmm_in[iter] <- result$ibs_est_train
        ibs_est_lcmm_out[iter] <- result$ibs_est_test
        ise_lcmm_in[iter] <- result$ISE_in
        ise_lcmm_out[iter] <- result$ISE_out
        CE_lcmm[iter] <- result$CE
      }
    }}
    
    # Save results for this parameter combination
    results_filename_nnem <- paste0("nnem_k3_results_n", n, "_c", gsub("\\.", "", c_str), ".RData")
    save(sigma_est, fixef_est, scale_est, shape_est, D_est, chol_est, omega_est, 
         pij, pij_true, CE, ibs_est_train, ibs_est_test, ISE_in, ISE_out, 
         MSE_y_out, MSE_y_in, file = results_filename_nnem)
    
    results_filename_lcmm <- paste0("jlcmm_k3_results_n", n, "_c", gsub("\\.", "", c_str), ".RData")
    save(sigma_est_lcmm, fixef_est_lcmm, scale_est_lcmm, shape_est_lcmm, D_est_lcmm, 
         chol_est_lcmm, omega_est_lcmm, pij_lcmm, pij_nnem, CE_lcmm, 
         ibs_est_lcmm_in, ibs_est_lcmm_out, ise_lcmm_in, ise_lcmm_out,
         MSE_y_out_lcmm, MSE_y_in_lcmm, file = results_filename_lcmm)
    
    cat("Completed n =", n, ", censoring rate =", c_rate, "\n")
  }
}

stopCluster(cl)