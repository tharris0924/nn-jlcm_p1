B <- 1000
n <- 1000
# theta_orig <- theta
sigma_est <- matrix(0, nrow = B)
fixef_est <- matrix(0, ncol = 6, nrow = B)
scale_est <- matrix(0, ncol = 3, nrow = B)
shape_est <- matrix(0, ncol = 3, nrow = B)
chol_est <- matrix(0, ncol = 3, nrow = B)
D_est <- matrix(0, ncol = 3, nrow = B)
omega_est <- matrix(0, ncol = 2, nrow = B)
CE <- ibs_est_lcmm_in <- ibs_est_lcmm_out <- ise_lcmm_in <- ise_lcmm_out <- numeric(
  B
)
library(matrixcalc)
options(width = 80)
pij <- pij_nnem <- c_true <- pred_c <- list()
pij_true <- list()
n <- 1000
for (i in seq_len(B)) {
  i <- 1
  wkdir <- paste0("NNEM/Example_k-3/k=3_Results_n=", n, "_c0.05")
  file <- paste0(wkdir, "/jlcmm_fit", i, n, ".RData")
  load(file)
  wkdir <- paste0("NNEM/Example_k-3/k=3_Results_n=", n, "_c0.05")
  file <- paste0(wkdir, "/nnem_size", i, n, ".RData")
  load(file)
  if (!is.null(result) && (!is.null(jlcmm_fit))) {
    if (!is.na(jlcmm_fit$best[1] )) {
      est <- jlcmm_fit$best
      obj <- jlcmm_fit
      nclasses <- ncol(obj$pprob) - 2
      coefs <- obj$best
      coefend <- min(which(grepl('Weibull', names(coefs)))) - 1
      coefs_multilogit <- matrix(coefs[1:coefend], nrow = nclasses - 1)
      tmpX <- cbind(1, result$x)
      linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
      exps = exp(cbind(linearval, 0))
      pi_logit <- exps / rowSums(exps)
      pij[i] <- list(pi_logit)
      # pij_true[i] <- list(result$theta_orig$pj)
      # set.seed(i)
      # C <- result$data_long$g
      # true_components <- numeric(500)

      true_components <- numeric(B)
      for (j in seq_len(n)) {
        red_c <- result$data_long[result$data_long$ID == j, ]
        true_components[j] <- unique(red_c$g)
      }
      true_mixing_prop <- function(x) {
        pi1 <- (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4)) * as.numeric(x < 3)
        pi2 <- (1 - (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4))) *
          as.numeric(x < 0)
        pi3 <- (2 * exp(-0.1 * (x - 6)^4)) /
          (1 + exp(-0.1 * (x - 6)^4)) *
          as.numeric(x >= 3)
        pi4 <- 1 - (pi1 + pi3)
        return(cbind(pi1, pi3, pi4))
      }
      # library(dplyr)
      pi_true <- true_mixing_prop(result$x)
      CE[i] <- mean(jlcmm_fit$pprob$class != true_components)
      pij_nnem[i] <- list(result$theta$pj)
      # pij_true[i] <- list(pi_true)
      # pij[i] <- list(theta$pj)
      pij_true[i] <- list(pi_true)
      sigma_est[i] <- est["stderr"]
      fixef_est[i, ] <- est[c(13, 16, 12, 15, 11, 14)]
      idxs <- which(grepl('Weibull', names(coefs)))
      shape_est[i, ] <- exp(est[c(10, 8, 6)])
      scale_est[i, ] <- exp(est[c(9, 7, 5)])^(-1 / exp(est[c(10, 8, 6)]))
      idxs <- which(grepl('varcov', names(est)))

      # omega_est[i,] <- matrix(sqrt((theta$omegas)^2), nrow = 1, byrow=T)
      chol_est[i, ] <- est[idxs]
      D_est[i, ] <- coefs[idxs]
      idxs <- which(grepl('varprop', names(est)))
      omega_est[i, ] <- est[idxs]
      ise_lcmm_in[i] <- get_lcmm_ISE(
        model = jlcmm_fit,
        data = result$data_wide,
        groups = result$x,
        theta_orig$shape,
        theta_orig$scale,
        pi_true
      )
      result$data_wide$x <- result$x
      ibs_est_lcmm_in[i] <- get_lcmm_brier_score(
        jlcmm_fit,
        result$data_wide,
        result$data_wide,
        NULL,
        ~x
      )

      #want to change the data
      wkdir <- paste0("NNEM/Example_k-3/k=3_Results_n=", n, "_c0.05")
      if (i == 1000) {
        file <- paste0(wkdir, "/nnem_size", 1, n, ".RData")
      } else {
        file <- paste0(wkdir, "/nnem_size", i + 1, n, ".RData")
      }
      load(file)
      if (!is.null(result)) {
        result$data_wide$x <- result$x
        ibs_est_lcmm_out[i] <- get_lcmm_brier_score(
          jlcmm_fit,
          result$data_wide,
          result$data_wide,
          NULL,
          ~x
        )
        # result$data_wide$x <- result$x
        ise_lcmm_out[i] <- get_lcmm_ISE(
          model = jlcmm_fit,
          data = result$data_wide,
          groups = result$x,
          theta_orig$shape,
          theta_orig$scale,
          pi_true
        )
      }
    }
    # theta <- result$theta
    # print("Error Variance:");theta$sigma_e;theta_orig$sigma_e;(bias_sigma <- abs(theta$sigma_e - theta_orig$sigma_e))
    extra <- nchar('||100%')
    width <- options()$width
    iter <- i
    maxit <- B
    step <- round(iter / maxit * (width - extra))
    text <- sprintf(
      '|%s%s|% 3s%%',
      strrep('=', step),
      strrep(' ', width - step - extra),
      round(iter / maxit * 100)
    )
    cat("\nSimulations Data Offload Progress:\n", text, '\n')
    # Sys.sleep(0.5)
    cat(if (iter == maxit) '\n' else '\014')
    # est_surv <- pred_surv_marg(shape_params=shape_est[i,], scale_params=shape_est[i,], mix.prob = pij[[i]], data_wide=result$data_wide)
    # times <- result$data_wide$event_time
    # deltas <- result$data_wide$delta
    # itx <- Surv(times, deltas)
    # ibs_est[i] <- sbrier(obj=itx,pred=est_surv)
  }
}
idxs_to_drop <- which(sigma_est > 0)
sigma_est <- sigma_est[idxs_to_drop]
scale_est <- scale_est[idxs_to_drop,]
shape_est <- shape_est[idxs_to_drop,]
fixef_est <- fixef_est[idxs_to_drop,]
D_est <- D_est[idxs_to_drop,]
chol_est <- chol_est[idxs_to_drop,]
omega_est <- omega_est[idxs_to_drop,]
# pij_bias <- pij_bias[idxs_to_drop,]
# mix.prop.bias <- colMeans(pij_bias)
lcmm_ce <- mean(CE, na.rm=T)

true_p <- MSE_p <- numeric(3)
for(k in seq_len(3)){
MAB_pi <- true_pi <- matrix(ncol=B, nrow = 1000)
for(i in seq_len(B)){
  # print(k)
if(!is.null(pij[i][[1]])){
MAB_pi[,i] <- pij[i][[1]][,k]
true_pi[,i] <- pij_true[i][[1]][,k]
}}
pi2 <- mean(abs(rowMeans(MAB_pi, na.rm = T) - rowMeans(true_pi, na.rm = T)));
pi3 <- mean(rowMeans(((MAB_pi) - (true_pi))^2, na.rm = T))
true_p[k] <- c(pi2)
MSE_p[k] <- pi3

}

# class_err <- mean(est_c)
mix.prop.bias <- true_p
# theta_orig <- result$theta_orig
theta_orig$betas <- c(t(theta_orig$betas))

ibs_est_lcmm_in_mean <- mean(ibs_est_lcmm_in[idxs_to_keep], na.rm = T)
ibs_est_lcmm_out_mean <- mean(ibs_est_lcmm_out[idxs_to_keep], na.rm = T)
ise_lcmm_in_mean <- mean(ise_lcmm_in, na.rm = T)
ise_lcmm_out_mean <- mean(ise_lcmm_out, na.rm = T)
print("Error Variance:");mean(sigma_est);theta_orig$sigma_e;#(bias_sigma <- abs(theta$sigma_e - theta_orig$sigma_e));mean(sigma_est)
print("Cholesky");theta_orig$cholesky; colMeans(chol_est); #(bias_cholesky <- abs(theta$cholesky - theta_orig$cholesky));colMeans(chol_est)
print("Omegas:");colMeans(sqrt(omega_est^2)); theta_orig$omegas;#(bias_omegas <- abs(sqrt((theta$omegas)^2) - theta_orig$omegas));colMeans(omega_est)
print("RE Variance-Covariance");colMeans(sqrt(omega_est^2))[1]^2 * colMeans(D_est);#vech(theta$D);(bias_D <- abs(vech(theta$D)-vech(theta_orig$D)));colMeans(D_est)
print("Longitudinal Fixef");colMeans(fixef_est); theta_orig$betas;#(bias_beta <- abs(theta$betas - theta_orig$betas));colMeans(fixef_est)
print("Weibull Scales");colMeans(scale_est); theta_orig$scale;#(bias_scale <- abs((theta$scale)-(theta_orig$scale)));colMeans(scale_est)
print("Weibull Shapes"); colMeans(shape_est); theta_orig$shape;# (bias_shape<- abs(theta$shape - theta_orig$shape));colMeans(shape_est)




bias_E <- abs(mean(sigma_est)-theta_orig$sigma_e);bias_E;
bias_scale <- abs(colMeans(scale_est)-matrix(theta_orig$scale, nrow = 1, byrow = T));bias_scale;colMeans(scale_est)
bias_shape <- abs(colMeans(shape_est)-matrix(theta_orig$shape, nrow = 1, byrow = T));bias_shape;colMeans(shape_est)
bias_fixef <- abs(colMeans(fixef_est)-matrix(theta_orig$betas, nrow = 1, byrow = T));bias_fixef; MAB_fixef <- (1/(3*(2+1)))*sum(bias_fixef); MAB_fixef;colMeans(fixef_est)
bias_D <- abs(colMeans(sqrt(omega_est^2))[1]^2 *colMeans(D_est)-matrix(vech(theta_orig$D), nrow = 1, byrow = T));bias_D;colMeans(sqrt(omega_est^2))[1]^2 * colMeans(D_est);
# bias_cholesky <- abs(colMeans(chol_est)-matrix(theta_orig$cholesky[lower.tri(theta$cholesky, diag = T)], nrow = 1, byrow = T));bias_cholesky;colMeans(chol_est)
bias_omegas <- abs(1/colMeans(sqrt(omega_est^2))-matrix(theta_orig$omegas, nrow = 1, byrow = T));bias_omegas;1/colMeans(sqrt(omega_est^2))


# install.packages("caret")
# MSE
B <- length(sigma_est)
matE <- rep(theta_orig$sigma_e,B)
MSE_sigma <- 1/B * t(sqrt(sigma_est)-matE)%*%(sqrt(sigma_est)-matE)


matB <- matrix(rep(theta_orig$betas,B), nrow=B, ncol=6, byrow = T)
MSE_beta <- 1/B * sum(colSums((fixef_est-matB)^2))

library(latex2exp)
library(kableExtra)

library(Metrics)

mse(fixef_est, matB)
MSE_sigma

# sigma_e
abs_bias_E <- abs(mean(sqrt(sigma_est))-theta_orig$sigma_e)
MSE_E <- Metrics::mse(matE,sqrt(sigma_est))
# MAB_E <- mean(abs(sigma_est-theta_orig$sigma_e))
# MAB_E <- Metrics::mdae(matE, sqrt(sigma_est))
# install.packages('yardstick')

# beta


matB <- matrix(rep(theta_orig$betas,B), nrow=B, ncol=6, byrow = T)
# correction factor
# MSE_beta <- colMeans((fixef_est-matB)^2); sqrt(MSE_beta)
MAB_beta <-   colMeans(abs(fixef_est-matB)); MAB_beta; 
MSE_beta <- sapply(1:6, function(i) Metrics::mse(matB[,i], fixef_est[,i]))

#Tau (Cholesky of D)
# matT <- matrix(rep(theta_orig$cholesky[lower.tri(theta_orig$cholesky, diag=T)],B), nrow=B, ncol=3, byrow = T)
# ave_T <- colMeans(chol_est)
# bias_cholesky <- abs(ave_T-matrix(theta_orig$cholesky[lower.tri(theta$cholesky, diag = T)], nrow = 1, byrow = T));bias_cholesky;
# MSE_T <- colMeans((chol_est-matT)^2); sqrt(MSE_T)
# MAB_T <-   colMeans(abs(chol_est-matT)); MAB_T

#D
bias_D <- abs(colMeans(sqrt(omega_est^2))[1]^2 * colMeans(D_est)-matrix(vech(theta_orig$D), nrow = 1, byrow = T));bias_D;
ave_D <- colMeans(sqrt(omega_est^2))[1]^2 *colMeans(D_est)
matD <- matrix(rep(vech(theta_orig$D),B), nrow = B, ncol=3, byrow = T)
MSE_D <- sapply(1:3, function(i) Metrics::mse( matD[,i], colMeans(sqrt(omega_est^2))[1]^2 * D_est[,i]))
# MAB_D <-  colMeans(abs(D_est-matD)); MAB_D

#kappa (shapes)
matK <- matrix(rep(theta_orig$shape,B), nrow=B, ncol=3, byrow = T)
print("Weibull Shapes"); colMeans(shape_est); theta_orig$shape;# (bias_shape<- abs(theta$shape - theta_orig$shape));colMeans(shape_est)
MSE_K <- sapply(1:3, function(i) Metrics::mse(matK[,i], shape_est[,i]))
MAB_K <-  colMeans(abs(shape_est-matK)); MAB_K

# lambda (scales)
print("Weibull Scales");colMeans(scale_est); theta_orig$scale;#(bias_scale <- abs((theta$scale)-(theta_orig$scale)));colMeans(scale_est)
matL <- matrix(rep(theta_orig$scale,B), nrow=B, ncol=3, byrow = T)
MSE_L <- sapply(1:3, function(i) Metrics::mse(matL[,i], scale_est[,i])); MSE_L
# MAB_L <-  colMeans(abs(scale_est-matL)); MAB_L
# View(scale_est)
# omegas
ave_omegas <- 1/colMeans(sqrt(omega_est^2))
# bias_omegas <- abs(colMeans(omega_est)-matrix(theta_orig$omegas, nrow = 1, byrow = T));bias_omegas;
matO <- matrix(rep(theta_orig$omegas,B), nrow=B, ncol=2, byrow = T)
MSE_O <- colMeans((omega_est-matO)^2); sqrt(MSE_O)
MAB_O <-  colMeans(abs(omega_est-matO)); MAB_O
MSE_O <- sapply(1:2, function(i) mse(matO[,i], 1/sqrt((omega_est[,i])^2)))


# After running your bias.R script, create a comprehensive results table

# After running your bias.R script, create a comprehensive results table

library(knitr)
library(kableExtra)
library(xtable)

# Create comprehensive results table
results_table <- data.frame(
  Parameter = c(
    "$\\pi_1(x)$", "$\\pi_2(x)$", "$\\pi_3(x)$", 

    # Error variance
    "$\\sigma_e$",
    
    # Fixed effects (6 parameters)
    "$\\beta_{11}$", "$\\beta_{12}$", "$\\beta_{21}$", "$\\beta_{22}$", "$\\beta_{31}$", "$\\beta_{32}$",
    
    # Cholesky elements (3 parameters)
    "$d_{11}$", "$d_{21}$", "$d_{22}$",
    
    # Weibull shapes (3 parameters)
    "$\\kappa_1$", "$\\kappa_2$", "$\\kappa_3$",
    
    # Weibull scales (3 parameters)  
    "$\\lambda_1$", "$\\lambda_2$", "$\\lambda_3$",
    
    # Omega parameters (2 parameters)
    "$\\omega_1$", "$\\omega_2$"
  ),
  
  True_Value = c(
    # True values - you'll need to replace these with actual theta_orig values
    numeric(3),
    theta_orig$sigma_e,
    theta_orig$betas,
    vech(theta_orig$D),
    theta_orig$shape,
    theta_orig$scale,
    theta_orig$omegas
  ),
  
  Estimate = c(
    # Estimated values
    numeric(3),
    mean(sigma_est),
    colMeans(fixef_est),
    colMeans(sqrt(omega_est^2))[1]^2* colMeans(D_est),
    colMeans(shape_est),
    colMeans(scale_est),
    ave_omegas
  ),
  
  Bias = c(
    # Bias calculations
    true_p,
    (bias_E),
    bias_fixef,
    bias_D,
    bias_shape,
    bias_scale,
    bias_omegas
  ),
  
  # MAB = c(
  #   # Absolute bias
  #   MAB_E,
  #   MAB_beta,
  #   MAB_D,
  #   MAB_K,
  #   MAB_L,
  #   MAB_O
  # ),
  
  MSE = c(
    # MSE values
    MSE_p,
    MSE_sigma,
    MSE_beta,
    MSE_D,
    MSE_K,
    MSE_L,
    MSE_O
  ),
  
  RMSE = c(
    # Root MSE
    sqrt(MSE_p),
    sqrt(MSE_sigma),
    sqrt(MSE_beta),
    sqrt(MSE_D),
    sqrt(MSE_K),
    sqrt(MSE_L),
    sqrt(MSE_O)
  )
)

# Round numeric columns for better presentation
numeric_cols <- c("True_Value", "Estimate", "Bias", "MSE", "RMSE")
results_table[numeric_cols] <- lapply(results_table[numeric_cols], function(x) round(x, 3))
results_table$True_Value <- sprintf("%.10g", results_table$True_Value)

# Method 1: Using xtable for LaTeX export
print("=== XTABLE OUTPUT ===")
xtable_output <- xtable(results_table, 
                       caption = "Bias and MSE Results for JLCMM Parameter Estimates under 3 components (B=1000 simulations)",
                       label = "tab:bias_mse_results",
                       digits = 3)

# Print LaTeX code
print(xtable_output, 
      type = "latex",
      include.rownames = FALSE,
      sanitize.text.function = function(x) x,  # Preserve LaTeX math notation
      caption.placement = "top",
      table.placement = "H")

# Method 2: Using kable with kableExtra for more formatting options
# print("=== KABLE OUTPUT ===")
# kable_output <- kable(results_table, 
#                      format = "latex", 
#                      booktabs = TRUE,
#                      escape = FALSE,  # Important for LaTeX math notation
#                      caption = "Bias and MSE Results for Parameter Estimates (B=500 simulations)",
#                      col.names = c("Parameter", "True Value", "Estimate", "Bias", "|Bias|", "MSE", "RMSE")) %>%
#   kable_styling(latex_options = c("striped", "hold_position", "scale_down")) %>%
#   pack_rows("Error Variance", 1, 1) %>%
#   pack_rows("Fixed Effects", 2, 7) %>%
#   pack_rows("Cholesky Elements", 8, B) %>%
#   pack_rows("Weibull Shapes", 11, 13) %>%
#   pack_rows("Weibull Scales", 14, 16) %>%
#   pack_rows("Omega Parameters", 17, 18)

# print(kable_output)

# # Method 3: Create a summary table by parameter type
# summary_table <- data.frame(
#   Parameter_Type = c("Error Variance", "Fixed Effects", "Cholesky Elements", "Weibull Shapes", "Weibull Scales", "Omega Parameters"),
#   Count = c(1, 6, 3, 3, 3, 2),
#   Mean_Abs_Bias = c(
#     MAB_E,
#     mean(MAB_beta),
#     mean(MAB_T),
#     mean(MAB_K),
#     mean(MAB_L),
#     mean(MAB_O)
#   ),
#   Mean_MSE = c(
#     mse_E,
#     mean(MSE_beta),
#     mean(MSE_T),
#     mean(MSE_K),
#     mean(MSE_L),
#     mean(MSE_O)
#   ),
#   Mean_RMSE = c(
#     sqrt(mse_E),
#     mean(sqrt(MSE_beta)),
#     mean(sqrt(MSE_T)),
#     mean(sqrt(MSE_K)),
#     mean(sqrt(MSE_L)),
#     mean(sqrt(MSE_O))
#   )
# )

# # Round summary table
# summary_table[c("Mean_Abs_Bias", "Mean_MSE", "Mean_RMSE")] <- 
#   lapply(summary_table[c("Mean_Abs_Bias", "Mean_MSE", "Mean_RMSE")], function(x) round(x, 4))

# print("=== SUMMARY TABLE ===")
# summary_xtable <- xtable(summary_table,
#                         caption = "Summary of Bias and MSE by Parameter Type",
#                         label = "tab:summary_results",
#                         digits = 4)

# print(summary_xtable,
#       type = "latex",
#       include.rownames = FALSE,
#       caption.placement = "top",
#       table.placement = "!htbp")

# # Save tables to files for easy LaTeX inclusion
# write.table(print(xtable_output, type = "latex", include.rownames = FALSE, 
#                  sanitize.text.function = function(x) x),
#            file = "bias_mse_table.tex", 
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Tables created successfully!\n")
# cat("Main results table saved to: bias_mse_table.tex\n")
  # cat("You can include it in your LaTeX document with: \\input{bias_mse_table.tex}\n")
GOF <- rbind(lcmm_ce, ibs_est_lcmm_in_mean,ibs_est_lcmm_out_mean, ise_lcmm_in_mean, ise_lcmm_out_mean)
xtable_output <- xtable(GOF, 
                       caption = paste0("Metrics of JLCMM for sample size ",n),
                       label = "tab:bias_mse_results",
                       digits = 6)

# Print LaTeX code
print(xtable_output, 
      type = "latex",
      include.rownames = TRUE,
      sanitize.text.function = function(x) x,  # Preserve LaTeX math notation
      caption.placement = "top",
      table.placement = "H")

obj <- jlcmm_fit
nclasses <- ncol(obj$pprob)-2
coefs <- obj$best
coefend <- min(which(grepl('Weibull',names(coefs))))-1
coefs_multilogit <- matrix(coefs[1:coefend],nrow=nclasses-1)
tmpX <- cbind(1,result$x)
linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
exps=exp(cbind(linearval,0));
pi_logit <- exps/rowSums(exps)
# lines(x, pi_logit[,1])
# lines(x, pi_logit[,2])
# lines(x,sdsds$theta_orig$pj[,1])
# lines(x,sdsds$theta_orig$pj[,2])
# test_class <- t(apply(linearval, 1,function(x){ }))
pis <- NULL
true_p <- MSE_p <- numeric(3)
pis_logit <- pis_nnem <- NULL
for(k in seq_len(3)){
MAB_pi <- true_pi <-nnems_pis <- matrix(ncol=B, nrow = 1000)
for(i in seq_len(B)){
  # print(k)
if(!is.null(pij[i][[1]])){
MAB_pi[,i] <- pij[i][[1]][,k]
true_pi[,i] <- pij_true[i][[1]][,k]
nnems_pis[,i] <- pij_nnem[i][[1]][,k]
}}
pis <- cbind(pis,rowMeans(true_pi, na.rm = T))
pis_logit <- cbind(pis_logit, rowMeans(MAB_pi, na.rm = T))
pis_nnem <- cbind(pis_nnem,  rowMeans(nnems_pis, na.rm = T))
pi2 <- mean(abs(rowMeans(MAB_pi, na.rm = T) - rowMeans(true_pi, na.rm = T)));
pi3 <- mean(rowMeans(((MAB_pi) - (true_pi))^2, na.rm = T))
true_p[k] <- c(pi2)
MSE_p[k] <- pi3

}
result$pj <- as.data.frame(result$pj)
df <- data.frame(pi_true, pi_logit, result$pj,result$x)
x <- as.matrix(x)


# cols <- list(unlist(colnames(df)))
# # cols <- rep(list('true one', "true two", "jlcmm one", "jlcmm two", "nnem one", "nnem two"))
# tmpdf <- cbind(df[,i], df$x,i)

i <- 1
tmpdf <- cbind(df[,i], df$result.x,i)
df_long <- data.frame(tmpdf)
for(i in 2:ncol(df)-1){
  tmpdf <- cbind(df[,i], df$result.x,i)
  df_long <- rbind(df_long,tmpdf)
}
# df_long <- df_long|>rename(`Props`=`V1`, `x`=`V2`)

colnames(df_long) <- c("Props","x","i")
df_long$i <- factor(df_long$i)
library(forcats)

# my_factor <- factor(c("A", "B", "C", "A", "B"))

# Recode with forcats (part of tidyverse)
df_long$i <- fct_recode(df_long$i,
                       "True p1" = "1",
                       "True p2" = "2",
                       "True p3" = "3",
                      "JLCMM p1"="4",
                      "JLCMM p2"="5",
                      "JLCMM p3"="6",
                      "NNEM p1"="7",
                      "NNEM p2"="8",
                      "NNEM p3"="9")
library(ggthemes)
library(ggplot2)
library(scales)
cols <- hue_pal()(3)
ggplot(data = df_long,aes(x=x,y=Props, col=factor(i))) + geom_line() + theme_minimal() + labs(title = "Mixing Proporitons Comparisons for Three Component Model",
       x = "x", 
       y = latex2exp::TeX("\\pi(x)")) +theme(legend.title=element_blank())+
      scale_color_manual(values = c(cols[1], cols[1], cols[1], cols[2], cols[2], cols[2], cols[3], cols[3], cols[3]))

df <- data.frame(true_mixing_prop(x), pi_logit, result$theta$pj,result$x)
x <- as.matrix(x)

# pis <- 
# cols <- list(unlist(colnames(df)))
# # cols <- rep(list('true one', "true two", "jlcmm one", "jlcmm two", "nnem one", "nnem two"))
# tmpdf <- cbind(df[,i], df$x,i)

i <- 1
tmpdf <- cbind(df[,i], df$result.x,i)
df_long <- data.frame(tmpdf)
for(i in 2:ncol(df)-1){
  tmpdf <- cbind(df[,i], df$result.x,i)
  df_long <- rbind(df_long,tmpdf)
}
# df_long <- df_long|>rename(`Props`=`V1`, `x`=`V2`)

colnames(df_long) <- c("Props","x","i")
df_long$i <- factor(df_long$i)
library(forcats)

# my_factor <- factor(c("A", "B", "C", "A", "B"))

# Recode with forcats (part of tidyverse)
df_long$i <- fct_recode(df_long$i,
                       "True p1" = "1",
                       "True p2" = "2",
                       "True p3" = "3",
                      "JLCMM p1"="4",
                      "JLCMM p2"="5",
                      "JLCMM p3"="6",
                      "NNEM p1"="7",
                      "NNEM p2"="8",
                      "NNEM p3"="9")
library(ggthemes)
library(scales)
cols <- hue_pal()(3)
ggplot(data = df_long,aes(x=x,y=Props, col=factor(i))) + geom_line() + theme_minimal() + labs(title = "Averaged Mixing Proporitons Comparisons",
       x = "x", 
       y = latex2exp::TeX("\\pi(x)")) +theme(legend.title=element_blank())+
      scale_color_manual(values = c(cols[1], cols[1], cols[1], cols[2], cols[2], cols[2], cols[3], cols[3], cols[3]))



