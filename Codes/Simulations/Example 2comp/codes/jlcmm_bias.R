B <- 1000
n <- 1000
# theta_orig <- theta
sigma_est <- matrix(0, nrow = B)
fixef_est <- matrix(0, ncol=4, nrow=B)
scale_est <- matrix(0, ncol=2, nrow=B)
shape_est <- matrix(0, ncol = 2, nrow = B)
# chol_est <- matrix(0, ncol = 2, nrow = B)
D_est <- matrix(0, ncol = 1, nrow = B)
# omega_est <- matrix(0, ncol = 2,nrow = B)
pij <- c_true <- pred_c <- list()
pij_true <- list()
library(matrixcalc)
pij_bias <- matrix(B, ncol=2, nrow = B)

CE <- ibs_est_lcmm_in<- ibs_est_lcmm_out <- ise_lcmm_in <- ise_lcmm_out<-MSE_y_out <- MSE_y_in numeric(B)
options(width=80)
library(matrixcalc)
library(survival)
library(prodlim)
library(data.table)
library(doParallel)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
 stopCluster(cl)
options(width=80)
offload <- foreach(
      iter = 1:10,
      .packages = c("matrixcalc", "prodlim", "data.table", 
    "survival", "lcmm"),
      .export = c(
        "pred_surv_marg", "sbrier", "get_nnem_brier_score", "estep",
        "marginal_density_yis", "density_weibull", "f_td", "sigma_est"
      ),.combine = c, .verbose = T
    ) %dopar% {
wkdir <- paste0("NNEM/Example2/k=2_Results_n=",n,"_c0.05")
file <- paste0(wkdir,"/jlcmm_fit", iter, n,".RData")
load(file)
wkdir <- paste0("NNEM/Example2/k=2_Results_n=",n,"_c0.25")
file <- paste0(wkdir,"/nnem_size", iter, n,".RData")
load(file)
  if(!is.null(jlcmm_fit)&!is.null(result)){
  est <- lcmm::estimates(jlcmm_fit)
  obj <- jlcmm_fit
nclasses <- ncol(obj$pprob)-2
coefs <- obj$best
coefend <- min(which(grepl('Weibull',names(coefs))))-1
coefs_multilogit <- matrix(coefs[1:coefend],nrow=nclasses-1)
tmpX <- cbind(1,result$x)
linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
exps=exp(cbind(linearval,0));
pi_logit <- exps/rowSums(exps)
pij[iter] <- list(pi_logit)
pij_true[iter] <- list(result$theta_orig$pj)
set.seed(iter)
# C <- result$data_long$g
true_components <- numeric(n)

for(j in seq_len(n)){
  red_c <- result$data_long[result$data_long$ID==j,]
  true_components[j] <- unique(red_c$g)
}
CE[i] <- mean(jlcmm_fit$pprob$class  != true_components)
  
  sigma_est[iter] <- est["stderr"]
  fixef_est[iter,] <- est[c(8,10,7,9)]
  shape_est[iter,] <- exp(est[c(6,4)])
  scale_est[iter,] <- exp(est[c(5,3)])^(-1/exp(est[c(6,4)]))
  # omega_est[i,] <- matrix(sqrt((theta$omegas)^2), nrow = 1, byrow=T)
  D_est[iter] <- est["cholesky 1"]
      result$data_wide$x <- result$x
    mod <- jlcmm_fit
    predy <- (mod$pred$pred_ss)
    
    # predy_test_max <- apply(array(seq_along(predclass_test_max)),1,function(x){predy_test_raw[x,predclass_test_max[x]]})
    # predy_test_avg <- rowSums(predy_test_raw * predclass_test)

    MSE_y_in <- Metrics::mse(result$data_long$y,predy)
  # ise_lcmm_in[i] <- get_lcmm_ISE(model=jlcmm_fit, data=result$data_wide,groups=result$x)
  ibs_est_lcmm_in[iter] <- get_nnem_brier_score(shape_params=shape_est[iter,], scale_params=scale_est[iter,], mix.prob =jlcmm_fit$pprob[,3:4] , test_data=result$data_wide)$IBS$ibs


    #want to change the data
    wkdir <- paste0("NNEM/Example2/k=2_Results_n=",n,"_c0.05")
    if(iter==1000){
      file <- paste0(wkdir,"/nnem_size", 1, n,".RData")}else{
      file <- paste0(wkdir,"/nnem_size", iter+1, n,".RData")
    }
  load(file)
    if(!is.null(result)){
  result$data_wide$x <- result$x
  ibs_est_lcmm_out[iter] <-  get_nnem_brier_score(shape_params=shape_est[iter,], scale_params=scale_est[iter,], mix.prob =jlcmm_fit$pprob[,3:4] , test_data=result$data_wide)
  result$data_long$x <- result$x[result$data_long$ID]
   
theta <- list(pj=pi_logit, scale=scale_est[iter,], shape=shape_est[iter,], 
  betas=t(matrix(fixef_est[i,],2,2, byrow = T)), ranef=jlcmm_fit$predRE$intercept, cholesky=D_est[iter],
  D=D_est[iter], sigma_e=sigma_est[iter])
  pred_mixing <- estep(pi_logit, theta)$P
   predy_test_raw <- rowSums(pi_logit[result$data_long$ID,]*lcmm::predictY(jlcmm_fit,newdata=result$data_long)$pred)
      # ise_lcmm_out[i] <- get_lcmm_ISE(model=jlcmm_fit, data=result$data_wide,groups=result$x)
    }

    MSE_y_out <- Metrics::mse(result$data_long$y,predy_test_raw)


  # est_surv <- pred_surv_marg(shape_params=shape_est[i,], scale_params=scale_est[i,], mix.prob = pij[[i]], data_wide=result$data_wide)
  # times <- result$data_wide$event_time
  # deltas <- result$data_wide$delta
  # itx <- Surv(times, deltas)
  # ibs_est[i] <- get_nnem_brier_score(shape_params=shape_est[i,], scale_params=scale_est[i,], mix.prob = pij[[i]], test_data=result$data_wide)
  # ibs_est_alt[i] <-  get_nnem_brier_score(shape_params=shape_est[i,], scale_params=scale_est[i,], mix.prob = pij[[i]], test_data=result$data_wide)

  }
  # theta <- result$theta
  # print("Error Variance:");theta$sigma_e;theta_orig$sigma_e;(bias_sigma <- abs(theta$sigma_e - theta_orig$sigma_e))
  # extra <- nchar('||100%')
  # width <- options()$width
  # iter <- i
  # maxit <- B
  # step <- round(iter / maxit * (width - extra))
  # text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
  #                 strrep(' ', width - step - extra), round(iter / maxit * 100))
  # cat("\nSimulations Data Offload Progress:\n",text, '\n')
  # # Sys.sleep(0.5)
  # cat(if (iter == maxit) '\n' else '\014')
                list(list(
      sigma_est =  sigma_est[iter],
      fixef_est= fixef_est[iter,], 
  scale_est = scale_est[iter,],
  scale_est = shape_est[iter,],
  # omega_est[iter,],
  D_est=D_est[iter],
  pij=pij[iter],
  MSE_y_in= MSE_y_in[iter],
  MSE_y_out= MSE_y_in[iter],
  ibs_est_train=ibs_est_train[iter],
  ibs_est_test =ibs_est_test[iter],
   CE=CE[iter]
      # oob_idx = oob_idx,
      # oob_pred = oob_class
    ))
    }

# C-index computations
# est_surv <- pred_surv_marg(shape_params=shape_est[1,], scale_params=scale_est[1,], mix.prob = pij[[1]], data_wide=result$data_wide)
# times <- result$data_wide$event_time
# deltas <- result$data_wide$delta
# # itx <- Surv(times, deltas)
# est_ibs <- sbrier(obj=itx,pred=est_surv)
# (Cindex <- nftbart::Cindex(est_surv, times, deltas))

# mix.prop.bias <- colMeans(pij_bias)
lcmm_ce <- mean(CE, na.rm = T)
true_p <- MSE_p <- numeric(2)
for(k in seq_len(2)){
MAB_pi <- true_pi <- matrix(ncol=B, nrow = n)
for(i in seq_len(B)){
  if(!is.null( pij[i][[1]])){
  # print(k)
MAB_pi[,i] <- pij[i][[1]][,k]
true_pi[,i] <- pij_true[i][[1]][,k]
}}
pi2 <- mean(abs(rowMeans(MAB_pi, na.rm = T) - rowMeans(true_pi, na.rm = T)));
pi3 <- mean(rowMeans(((MAB_pi) - (true_pi))^2, na.rm = T))
true_p[k] <- c(pi2)
MSE_p[k] <- pi3

}


MSE_p
(mix.prop.bias <- true_p)
theta_orig <- result$theta_orig
theta_orig$betas <- c(t(theta_orig$betas))
idxs_to_keep <- which(sigma_est>0)
sigma_est <- sigma_est[idxs_to_keep]
D_est <- D_est[idxs_to_keep]
fixef_est <- fixef_est[idxs_to_keep,]
scale_est <- scale_est[idxs_to_keep,]
shape_est <- shape_est[idxs_to_keep,]
ibs_est_lcmm_in_mean <- mean(ibs_est_lcmm_in[idxs_to_keep], na.rm = T)
ibs_est_lcmm_out_mean <- mean(ibs_est_lcmm_out[idxs_to_keep], na.rm = T)
ise_lcmm_in_mean <- mean(ise_lcmm_in, na.rm = T)
ise_lcmm_out_mean <- mean(ise_lcmm_out, na.rm = T)

GOF <- rbind(lcmm_ce, ibs_est_lcmm_in_mean,ibs_est_lcmm_out_mean, ise_lcmm_in_mean, ise_lcmm_out_mean)
cat("\nError Variance:",mean(sigma_est),theta_orig$sigma_e)#(bias_sigma <- abs(theta$sigma_e - theta_orig$sigma_e));mean(sigma_est)
# print("Cholesky");theta_orig$cholesky; colMeans(chol_est); #(bias_cholesky <- abs(theta$cholesky - theta_orig$cholesky));colMeans(chol_est)
cat("\nRE Variance-Covariance",mean(D_est), sqrt(theta_orig$D[1]));#(bias_D <- abs(vech(theta$D)-vech(theta_orig$D)));colMeans(D_est)
# print("Omegas:");colMeans(sqrt(omega_est^2)); theta_orig$omegas;#(bias_omegas <- abs(sqrt((theta$omegas)^2) - theta_orig$omegas));colMeans(omega_est)
cat("\nLongitudinal Fixef",colMeans(fixef_est), theta_orig$betas);#(bias_beta <- abs(theta$betas - theta_orig$betas));colMeans(fixef_est)
cat("\nWeibull Scales",colMeans(scale_est), theta_orig$scale);#(bias_scale <- abs((theta$scale)-(theta_orig$scale)));colMeans(scale_est)
cat("\nWeibull Shapes", colMeans(shape_est), theta_orig$shape)# (bias_shape<- abs(theta$shape - theta_orig$shape));colMeans(shape_est)




bias_sigma <- abs(mean(sqrt(sigma_est))-theta_orig$sigma_e);bias_sigma;
bias_scale <- abs(colMeans(scale_est)-matrix(theta_orig$scale, nrow = 1, byrow = T));bias_scale;colMeans(scale_est)
bias_shape <- abs(colMeans(shape_est)-matrix(theta_orig$shape, nrow = 1, byrow = T));bias_shape;colMeans(shape_est)
bias_fixef <- abs(colMeans(fixef_est)-matrix(theta_orig$betas, nrow = 1, byrow = T));bias_fixef; MAB_fixef <- (1/(3*(2+1)))*sum(bias_fixef); MAB_fixef;colMeans(fixef_est)
bias_D <- abs(mean(D_est)-matrix(theta_orig$D, nrow = 1, byrow = T));bias_D;mean(D_est)
# bias_cholesky <- abs(colMeans(chol_est)-matrix(theta_orig$cholesky[lower.tri(theta$cholesky, diag = T)], nrow = 1, byrow = T));bias_cholesky;colMeans(chol_est)
# bias_omegas <- abs(colMeans(omega_est)-matrix(theta_orig$omegas, nrow = 1, byrow = T));bias_omegas;colMeans(omega_est)

library(lcmm)
# result$data_long$x <- result$x[result$data_long$ID] 

#     MSE_y_test_max <- mean(((predy_test_max) - data_test$y)^2)
#     MSE_y_test_avg <- mean((predy_test_avg - data_test$y)^2)
# install.packages("caret")
# MSE
B <- length(sigma_est)

matB <- matrix(rep(theta_orig$betas,B), nrow=B, ncol=4, byrow = T)
MSE_beta <- 1/B * sum(colSums((fixef_est-matB)^2))

library(latex2exp)
library(kableExtra)

library(Metrics)

mse(fixef_est, matB)
# MSE_sigma

# sigma_e
matE <- rep(theta_orig$sigma_e,B)
abs_bias_E <- abs(mean((sigma_est))-theta_orig$sigma_e)
MSE_E <- Metrics::mse(matE,(sigma_est))
# MAB_E <- mean(abs(sigma_est-theta_orig$sigma_e))
# MAB_E <- Metrics::mdae(matE, sqrt(sigma_est))
# install.packages('yardstick')

# beta


matB <- matrix(rep(theta_orig$betas,B), nrow=B, ncol=4, byrow = T)
# correction factor
MSE_beta <- sapply(1:4, function(i) Metrics::mse(matB[,i], fixef_est[,i]))
MAB_beta <-   colMeans(abs(fixef_est-matB)); MAB_beta; 
# MAB_beta <- sapply(1:4, function(i) Metrics::mdae(matB[,i], fixef_est[,i]))

#Tau (Cholesky of D)
# matT <- matrix(rep(theta_orig$cholesky[lower.tri(theta_orig$cholesky, diag=T)],B), nrow=B, ncol=3, byrow = T)
# ave_T <- colMeans(chol_est)
# bias_cholesky <- abs(ave_T-matrix(theta_orig$cholesky[lower.tri(theta$cholesky, diag = T)], nrow = 1, byrow = T));bias_cholesky;
# MSE_T <- colMeans((chol_est-matT)^2); sqrt(MSE_T)
# MAB_T <-   colMeans(abs(chol_est-matT)); MAB_T

#D
bias_D <- abs(mean((D_est))-matrix(sqrt(theta_orig$D[1]), nrow = 1, byrow = T));bias_D;
ave_D <- mean((D_est))
matD <- matrix(theta_orig$D, nrow = B, ncol=1, byrow = T)
# MSE_D <- mean((sqrt(D_est)-sqrt(matD))^2); (MSE_D)
# (MSE_D <- var(sqrt(D_est))+(bias_D)^2)
MSE_D <- sapply(1:1, function(i) Metrics::mse(sqrt(matD), (D_est)))
# MAB_D <-  mean(abs(sqrt(D_est)-sqrt(matD))); MAB_D

#kappa (shapes)
matK <- matrix(rep(theta_orig$shape,B), nrow=B, ncol=2, byrow = T)
print("Weibull Shapes"); colMeans(shape_est); theta_orig$shape;# (bias_shape<- abs(theta$shape - theta_orig$shape));colMeans(shape_est)
# MSE_K <- colMeans((shape_est-matK)^2);
MSE_K <- sapply(1:2, function(i) Metrics::mse(matK[,i], shape_est[,i]))

MAB_K <-  colMeans(abs(shape_est-matK)); MAB_K

# lambda (scales)
print("Weibull Scales");colMeans(scale_est); theta_orig$scale;#(bias_scale <- abs((theta$scale)-(theta_orig$scale)));colMeans(scale_est)
matL <- matrix(rep(theta_orig$scale,B), nrow=B, ncol=2, byrow = T)
MSE_L <- colMeans((scale_est-matL)^2); sqrt(MSE_L)
MAB_L <-  colMeans(abs(scale_est-matL)); MAB_L
(MSE_L <-sapply(1:2, function(i) var(scale_est[,i])+(bias_scale[,i])^2))
MSE_L <- sapply(1:2, function(i) Metrics::mse(matL[,i], scale_est[,i]))
MSE_L
# View(scale_est)
# View(scale_est)
# omegas
# ave_omegas <- colMeans(omega_est)
# bias_omegas <- abs(colMeans(omega_est)-matrix(theta_orig$omegas, nrow = 1, byrow = T));bias_omegas;
# matO <- matrix(rep(theta_orig$omegas,B), nrow=B, ncol=2, byrow = T)
# MSE_O <- colMeans((omega_est-matO)^2); sqrt(MSE_O)
# MAB_O <-  colMeans(abs(omega_est-matO)); MAB_O

# After running your bias.R script, create a comprehensive results table

# After running your bias.R script, create a comprehensive results table

library(knitr)
library(kableExtra)
library(xtable)

# Create comprehensive results table
results_table <- data.frame(
  Parameter = c(

    # Mixing proportions 
    "$\\pi_1(x)$", "$\\pi_2(x)$",
    
    # Error variance
    "$\\sigma_e$",
    
    # Fixed effects (6 parameters)
    "$\\beta_{11}$", "$\\beta_{12}$", "$\\beta_{21}$", "$\\beta_{22}$",
    
    # Cholesky elements (3 parameters)
    "$d_{11}$",
    
    # Weibull shapes (3 parameters)
    "$\\kappa_1$", "$\\kappa_2$",
    
    # Weibull scales (3 parameters)  
    "$\\lambda_1$", "$\\lambda_2$"
    
  ),
  
  True_Value = c(
    # True values - you'll need to replace these with actual theta_orig values
    numeric(2),
    theta_orig$sigma_e,
    theta_orig$betas,
    sqrt(theta_orig$D[1]),
    theta_orig$shape,
    theta_orig$scale
    # theta_orig$omegas
  ),
  
  Estimate = c(
    # Estimated values
    numeric(2),
    mean((sigma_est)),
    colMeans(fixef_est),
    mean((D_est)),
    colMeans(shape_est),
    colMeans(scale_est)
    # colMeans(omega_est)
  ),
  
  Bias = c(
    # Bias calculations
    mix.prop.bias,
    abs_bias_E,
    bias_fixef,
    bias_D[1],
    bias_shape,
    bias_scale
  ),
  
  # MAB = c(
  #   # Absolute bias
  #   MAB_E,
  #   MAB_beta,
  #   MAB_D,
  #   MAB_K,
  #   MAB_L
  # ),
  
  MSE = c(
    # MSE values
    MSE_p,
    MSE_E,
    MSE_beta,
    MSE_D,
    MSE_K,
    MSE_L
  ),
  
  RMSE = c(
    # Root MSE
    sqrt(MSE_p),
    sqrt(MSE_E),
    sqrt(MSE_beta),
    sqrt(MSE_D),
    sqrt(MSE_K),
    sqrt(MSE_L)
  )
)

# Round numeric columns for better presentation
numeric_cols <- c("True_Value", "Estimate", "Bias", "MSE", "RMSE")
results_table[numeric_cols] <- lapply(results_table[numeric_cols], function(x) round(x, 3))
results_table$True_Value <- sprintf("%.10g", results_table$True_Value)

# Method 1: Using xtable for LaTeX export
library(xtable)
print("=== XTABLE OUTPUT ===")
xtable_output <- xtable(results_table, 
                       caption = "Bias and MSE Results for JLCMM Parameter Estimates (B=500 simulations)",
                       label = "tab:bias_mse_results",
                       digits = 3)

# Print LaTeX code
print(xtable_output, 
      type = "latex",
      include.rownames = FALSE,
      sanitize.text.function = function(x) x,  # Preserve LaTeX math notation
      caption.placement = "top",
      table.placement = "H")


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