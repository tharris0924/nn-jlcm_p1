# Define parameter combinations
sample_sizes <- c(100, 500, 1000)[1]
censoring_rates <- c(0.05, 0.25, 0.5)

library(doParallel)
library(lcmm)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
n <- 100
# Outer loop for sample sizes
for(n in sample_sizes) {
  # Inner loop for censoring rates
  for(c_rate in censoring_rates) {
    
    cat("Processing n =", n, ", censoring rate =", c_rate, "\n")
    
    # Initialize result storage for this combination
    B <- 1000
    
    # NNEM results
    sigma_est <- matrix(0, nrow = B)
    fixef_est <- matrix(0, ncol=4, nrow=B)
    scale_est <- matrix(0, ncol=2, nrow=B)
    shape_est <- matrix(0, ncol = 2, nrow = B)
    D_est <- matrix(0, ncol = 1, nrow = B)
    pij <- c_true <- pred_c <- list()
    pij_true <- list()
    CE <- ibs_est_train <- ibs_est_test <- ISE_in <- ISE_out <- MSE_y_out <- MSE_y_in <- numeric(B)
    
    # JLCMM results
    sigma_est_lcmm <- matrix(0, nrow = B)
    fixef_est_lcmm <- matrix(0, ncol=4, nrow=B)
    scale_est_lcmm <- matrix(0, ncol=2, nrow=B)
    shape_est_lcmm <- matrix(0, ncol = 2, nrow = B)
    D_est_lcmm <- matrix(0, ncol = 1, nrow = B)
    pij_lcmm <- list()
    CE_lcmm <- ibs_est_lcmm_in <- ibs_est_lcmm_out <- MSE_y_out_lcmm <- MSE_y_in_lcmm <- numeric(B)
    
    # Convert censoring rate to string format for file paths
    c_str <- if(c_rate == 0.05) "0.05" else if(c_rate == 0.25) "0.25" else "0.5"
    
    # Your parallel processing
    offload <- foreach(
      iter = 1:1000,
      .packages = c("matrixcalc", "prodlim", "data.table", "survival", "SurvMetrics", "lcmm"),
      .export = c(
        "pred_surv_marg", "sbrier", "get_nnem_brier_score", "estep",
        "marginal_density_yis", "density_weibull", "f_td", "n", "c_str"
      ),
      .combine = c, 
      .verbose = TRUE, .errorhandling = "pass"
    ) %dopar% {
      
      # Load JLCMM results
      wkdir_jlcmm <- paste0("NNEM/Example2/k=2_Results_n=", n, "_c", c_str)
      file_jlcmm <- paste0(wkdir_jlcmm, "/jlcmm_fit", iter, n, ".RData")
      
      # Load NNEM results  
      wkdir_nnem <- paste0("NNEM/Example2/k=2_Results_n=", n, "_c", c_str)
      file_nnem <- paste0(wkdir_nnem, "/nnem_size", iter, n, ".RData")
      
      # Initialize result containers
      nnem_result <- NULL
      jlcmm_result <- NULL
      
      # Process NNEM results
      if(file.exists(file_nnem)) {
        load(file_nnem)
        theta <- result$theta
        
        if(!is.null(theta)){
          # Full NNEM computation code
          local_sigma_est <- theta$sigma_e^2
          local_fixef_est <- matrix(theta$betas, nrow=1, byrow=T)
          local_scale_est <- matrix(theta$scale, nrow = 1, byrow = T) 
          local_shape_est <- matrix(theta$shape, nrow=1, byrow = T)
          local_D_est <- theta$D
          local_pij <- theta$pj
          local_pij_true <- result$theta_orig$pj
          
          true_components <- numeric(n)
          for(j in seq_len(n)){
            red_c <- result$data_long[result$data_long$ID==j,]
            true_components[j] <- unique(red_c$g)
          }
          
          # Prediction code for NNEM
          pred_y <- matrix(0, nrow(result$data_long), 2)
          Z <- model.matrix(~1, data = result$data_long)
          X <- model.matrix(~measure_time, data=result$data_long)
          
          data <- result$data_long
          data_wide <- result$data_wide
          
          f_td <- function(shape, scale, t, d) {
            ((dweibull(t, shape=shape, scale = scale)^d)) * 
            (pweibull(t, shape = shape, scale = scale, lower.tail = F)^((1-d))) 
          }
          # theta$betas <- t(theta$betas)
          comps <- estep(prior_mixing = theta$pj, theta)
          comp.p <- comps$P
          n_timepoints <- table(result$data_long$ID)     
          bi <- list()
          
          m <- length(theta$shape)
          q <- nrow(theta$D)
          
          for(j in 1:m){
            sigma_g <- result$theta$sigma_e
            D_g <- theta$D
            beta_g <- theta$betas[,j]
            
            bigs <- matrix(nrow = n, ncol = q)
            for(k in 1:n){
              n_i <- n_timepoints[k]
              sim <- data
              subset <- sim[sim$ID==k,]
              Z_i <- Z[sim$ID==k, , drop = FALSE]
              X_i <- X[sim$ID==k, ,drop=FALSE]
              p <- ncol(X_i)
              Xi_beta <- X_i %*% beta_g
              q <- ncol(Z_i)
              Y_i <- subset$y
              sigma <- sigma_g
              
              Vi <- as.matrix(Z_i%*% D_g%*% t(Z_i) + sigma_g^2 * diag(n_i))
              big <- D_g %*% t(Z_i) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
              bigs[k] <- big
            }
            bi[j] <- list(bigs)
          }
          
          for(g in 1:2){
            mu <- NULL
            beta_tmp <- theta$betas[,g]
            ranef_tmp <- bi[[g]]
            
            for (k in 1:n) {
              idx <- which(result$data_long$ID == k)
              Z_i <- Z[idx, , drop = FALSE]
              X_i <- X[idx, ,drop=FALSE]
              XBg <- (X_i %*% beta_tmp + Z_i%*%ranef_tmp[k])
              mu <- c(mu, XBg)
            }
            pred_y[,g] <- mu
          }
          
          preds <- rowSums(pred_y * comp.p[result$data_long$ID,] )
          
          local_MSE_y_in <- Metrics::mse(actual = (result$data_long$y), predicted = (preds))
          local_CE <- mean(result$C != true_components)
          
          # data <- result$data_long
          # data_wide <- result$data_wide
          # result$theta_orig$betas <- t(result$theta_orig$betas)
          # true_comp <- estep(result$theta_orig$pj, result$theta_orig)$P
          
          local_ibs_est_train <- get_nnem_brier_score(
            shape_params=local_shape_est, 
            scale_params=local_scale_est, 
            mix.prob = comp.p, 
            test_data=result$data_wide
          )$IBS$ibs
          
          # Out-of-sample evaluation
          if(iter==1000){
            file_test <- paste0(wkdir_nnem,"/nnem_size", 1, n,".RData")
          } else {
            file_test <- paste0(wkdir_nnem,"/nnem_size", iter+1, n,".RData")
          }
          
          local_MSE_y_out <- NA
          local_ibs_est_test <- NA
          
          if(file.exists(file_test)) {
            load(file_test)
            if(!is.null(result)){
              data <- result$data_long
              data_wide <- result$data_wide
              result$theta$betas <- t(result$theta$betas)
              true_comp <- estep(result$theta$pj, result$theta)$P
              
              # Repeat prediction process for out-of-sample data
              pred_y <- matrix(0, nrow(result$data_long), 2)
              Z <- model.matrix(~1, data = result$data_long)
              X <- model.matrix(~measure_time, data=result$data_long)
              
              n_timepoints <- table(result$data_long$ID)     
              bi <- list()
              
              for(j in 1:m){
                sigma_g <- theta$sigma_e # Use theta from training data
                D_g <- theta$D
                beta_g <- theta$betas[,j]
                
                bigs <- matrix(nrow = n, ncol = q)
                for(k in 1:n){
                  n_i <- n_timepoints[k]
                  sim <- data
                  subset <- sim[sim$ID==k,]
                  Z_i <- Z[sim$ID==k, , drop = FALSE]
                  X_i <- X[sim$ID==k, ,drop=FALSE]
                  p <- ncol(X_i)
                  Xi_beta <- X_i %*% beta_g
                  q <- ncol(Z_i)
                  Y_i <- subset$y
                  sigma <- sigma_g
                  
                  Vi <- as.matrix(Z_i%*% D_g%*% t(Z_i) + sigma_g^2 * diag(n_i))
                  big <- D_g %*% t(Z_i) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
                  bigs[k] <- big
                }
                bi[j] <- list(bigs)
              }
              
              for(g in 1:2){
                mu <- NULL
                beta_tmp <- theta$betas[,g]
                ranef_tmp <- bi[[g]]
                
                for (k in 1:n) {
                  idx <- which(result$data_long$ID == k)
                  Z_i <- Z[idx, , drop = FALSE]
                  X_i <- X[idx, ,drop=FALSE]
                  XBg <- (X_i %*% beta_tmp + Z_i%*%ranef_tmp[k])
                  mu <- c(mu, XBg)
                }
                pred_y[,g] <- mu
              }
              
              preds <- rowSums(pred_y * true_comp[result$data_long$ID,] )
              
              local_MSE_y_out <- Metrics::mse(actual = (result$data_long$y), predicted = (preds))
              local_ibs_est_test <- get_nnem_brier_score(
                shape_params=local_shape_est, 
                scale_params=local_scale_est, 
                mix.prob = true_comp, 
                test_data=result$data_wide
              )$IBS$ibs
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
            pij = local_pij,
            pij_true = local_pij_true,
            MSE_y_in = local_MSE_y_in,
            MSE_y_out = local_MSE_y_out,
            ibs_est_train = local_ibs_est_train,
            ibs_est_test = local_ibs_est_test,
            CE = local_CE
          )
        }
      }
      
      # Process JLCMM results
      if(file.exists(file_jlcmm)) {
        jlcmm_result <- tryCatch({
        load(file_jlcmm)  # This loads jlcmm_fit
        load(file_nnem)   # Also need result object for comparison
        
        if(!is.null(jlcmm_fit) & !is.null(result)){
          est <- jlcmm_fit$best
          obj <- jlcmm_fit
          nclasses <- ncol(obj$pprob)-2
          coefs <- obj$best
          coefend <- min(which(grepl('Weibull',names(coefs))))-1
          coefs_multilogit <- matrix(coefs[1:coefend], nrow=nclasses-1)
          tmpX <- cbind(1, result$x)
          linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
          exps <- exp(cbind(linearval, 0))
          pi_logit <- exps/rowSums(exps)
          
          # True components
          true_components <- numeric(n)
          for(j in seq_len(n)){
            red_c <- result$data_long[result$data_long$ID==j,]
            true_components[j] <- unique(red_c$g)
          }
          
          local_CE_lcmm <- mean(jlcmm_fit$pprob$class != true_components)
          local_sigma_est_lcmm <- est["stderr"]
          local_fixef_est_lcmm <- est[c(8,10,7,9)]
          local_shape_est_lcmm <- exp(est[c(6,4)])
          local_scale_est_lcmm <- exp(est[c(5,3)])^(-1/exp(est[c(6,4)]))
          local_D_est_lcmm <- est["varcov 1"]
          
          # In-sample evaluation
          result$data_wide$x <- result$x
          mod <- jlcmm_fit
          predy <- mod$pred$pred_ss
          local_MSE_y_in_lcmm <- Metrics::mse(result$data_long$y, predy)
          local_ibs_est_lcmm_in <- get_nnem_brier_score(
            shape_params=local_shape_est_lcmm, 
            scale_params=local_scale_est_lcmm, 
            mix.prob=jlcmm_fit$pprob[,3:4], 
            test_data=result$data_wide
          )$IBS$ibs
          
          # Out-of-sample evaluation
          if(iter==1000){
            file_test <- paste0(wkdir_nnem,"/nnem_size", 1, n,".RData")
          } else {
            file_test <- paste0(wkdir_nnem,"/nnem_size", iter+1, n,".RData")
          }
          
          local_MSE_y_out_lcmm <- 0
          local_ibs_est_lcmm_out <- 0
          
          if(file.exists(file_test)) {
            load(file_test)
            if(!is.null(result)){
              result$data_wide$x <- result$x
              result$data_long$x <- result$x[result$data_long$ID]
              tmpX <- cbind(1, result$x)
              linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
              exps <- exp(cbind(linearval, 0))
              pi_logit <- exps/rowSums(exps)
              theta <- list(
                pj=pi_logit, 
                scale=local_scale_est_lcmm, 
                shape=local_shape_est_lcmm,
                betas=t(matrix(local_fixef_est_lcmm, 2, 2, byrow = T)), 
                ranef=jlcmm_fit$predRE$intercept, 
                cholesky=local_D_est_lcmm,
                D=local_D_est_lcmm, 
                sigma_e=local_sigma_est_lcmm
              )
                        f_td <- function(shape, scale, t, d) {
            ((dweibull(t, shape=shape, scale = scale)^d)) * 
            (pweibull(t, shape = shape, scale = scale, lower.tail = F)^((1-d))) 
          }
              data <- result$data_long
              data_wide <- result$data_wide
              
      
              pred_mixing <- estep(pi_logit, theta)$P
              local_ibs_est_lcmm_out <- get_nnem_brier_score(
                shape_params=local_shape_est_lcmm, 
                scale_params=local_scale_est_lcmm, 
                mix.prob=pred_mixing, 
                test_data=result$data_wide
              )$IBS$ibs
              

              predy_test_raw <- rowSums(pred_mixing[result$data_long$ID,] * 
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
            pij = pi_logit,
            MSE_y_in = local_MSE_y_in_lcmm,
            MSE_y_out = local_MSE_y_out_lcmm,
            ibs_est_train = local_ibs_est_lcmm_in,
            ibs_est_test = local_ibs_est_lcmm_out,
            CE = local_CE_lcmm
          )
        }
          }, error=function(e){
          warning(paste("JLCMM computation failed for iter =", iter, ". Error:", e$message))
          return(NULL)
          })
        }
      
      # Return both results
      return(list(nnem = nnem_result, jlcmm = jlcmm_result))
    }
    
    # Process results and assign to appropriate variables
    for(result_pair in offload) {
      # Process NNEM results
      if(!is.null(result_pair) && result_pair$method=="NNEM") {
        result <- result_pair
        iter <- result$iter
        sigma_est[iter] <- result$sigma_est
        fixef_est[iter,] <- result$fixef_est
        scale_est[iter,] <- result$scale_est
        shape_est[iter,] <- result$shape_est
        D_est[iter] <- result$D_est
        pij[iter] <- list(result$pij)
        pij_true[iter] <- list(result$pij_true)
        MSE_y_in[iter] <- result$MSE_y_in
        MSE_y_out[iter] <- result$MSE_y_out
        ibs_est_train[iter] <- result$ibs_est_train
        ibs_est_test[iter] <- result$ibs_est_test
        CE[iter] <- result$CE
      }
      
      # Process JLCMM results
      if(!is.null(result_pair) && result_pair$method=="JLCMM") {
        result <- result_pair 
        iter <- result$iter
        sigma_est_lcmm[iter] <- result$sigma_est
        fixef_est_lcmm[iter,] <- result$fixef_est
        scale_est_lcmm[iter,] <- result$scale_est
        shape_est_lcmm[iter,] <- result$shape_est
        D_est_lcmm[iter] <- result$D_est
        pij_lcmm[iter] <- list(result$pij)
        MSE_y_in_lcmm[iter] <- result$MSE_y_in
        MSE_y_out_lcmm[iter] <- result$MSE_y_out
        ibs_est_lcmm_in[iter] <- result$ibs_est_train
        ibs_est_lcmm_out[iter] <- result$ibs_est_test
        CE_lcmm[iter] <- result$CE
      }
    }
    
    # Save results for this parameter combination
    results_filename_nnem <- paste0("nnem_results_n", n, "_c", gsub("\\.", "", c_str), ".RData")
    save(sigma_est, fixef_est, scale_est, shape_est, D_est, pij, pij_true,
         CE, ibs_est_train, ibs_est_test, ISE_in, ISE_out, MSE_y_out, MSE_y_in,
         file = results_filename_nnem)
    
    results_filename_lcmm <- paste0("jlcmm_results_n", n, "_c", gsub("\\.", "", c_str), ".RData")
    save(sigma_est_lcmm, fixef_est_lcmm, scale_est_lcmm, shape_est_lcmm, D_est_lcmm, pij_lcmm,
         CE_lcmm, ibs_est_lcmm_in, ibs_est_lcmm_out, MSE_y_out_lcmm, MSE_y_in_lcmm,
         file = results_filename_lcmm)
    
    cat("Completed n =", n, ", censoring rate =", c_rate, "\n")
  }
}

stopCluster(cl)