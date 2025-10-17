B <- 1000
n <- 1000
library(SurvMetrics)
library('survival')
# theta_orig <- theta
sigma_est <- matrix(0, nrow = B)
fixef_est <- matrix(0, ncol=4, nrow=B)
scale_est <- matrix(0, ncol=2, nrow=B)
shape_est <- matrix(0, ncol = 2, nrow = B)
# chol_est <- matrix(0, ncol = 2, nrow = B)
D_est <- matrix(0, ncol = 1, nrow = B)
# omega_est <- matrix(0, ncol = 2,nrow = B)
source("C:/Users/th240/Documents/Ph.D. Statistics/Codes/NNEM/Example2/em_utils copy.R")
source("c:/Users/th240/Documents/Ph.D. Statistics/Codes/NNEM/Example2/estep.R")
pij <- c_true <- pred_c <- list()
pij_true <- list()
CE <- ibs_est_train <- ibs_est_test <-ISE_in <- ISE_out <- MSE_y_out <- MSE_y_in <- numeric(B)
library(matrixcalc)
library(survival)
library(prodlim)
library(data.table)
library(doParallel)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
# stopCluster(cl)
options(width=80)
# Define batch size
# batch_size <- 100
# Define parameter combinations
sample_sizes <- c(100, 500, 1000)
censoring_rates <- c(0.05, 0.25, 0.5)
n_total <- 1000
# batch_size <- 50
# batch_indices <- split(1:n_total, ceiling((1:n_total) / batch_size))
# Create batches
# print("First batch contents:")
# print(batch_indices[[1]])
# Process batches
library(tictoc)
tic()
for(n in sample_sizes) {
  # Inner loop for censoring rates
  for(c_rate in censoring_rates) {
    B <- 1000
    c_rate <- 0.05
    n <- 100
library(SurvMetrics)
library('survival')
# theta_orig <- theta
sigma_est <- matrix(0, nrow = B)
fixef_est <- matrix(0, ncol=4, nrow=B)
scale_est <- matrix(0, ncol=2, nrow=B)
shape_est <- matrix(0, ncol = 2, nrow = B)
# chol_est <- matrix(0, ncol = 2, nrow = B)
D_est <- matrix(0, ncol = 1, nrow = B)
# omega_est <- matrix(0, ncol = 2,nrow = B)
source("C:/Users/th240/Documents/Ph.D. Statistics/Codes/NNEM/Example2/em_utils copy.R")
source("c:/Users/th240/Documents/Ph.D. Statistics/Codes/NNEM/Example2/estep.R")
pij <- c_true <- pred_c <- list()
pij_true <- list()
CE <- ibs_est_train <- ibs_est_test <-ISE_in <- ISE_out <- MSE_y_out <- MSE_y_in <- numeric(B)
offload <- foreach(
  iter=1:10,
  .packages = c("matrixcalc", "prodlim", "data.table", "survival"),
  .export = c(
    "pred_surv_marg", "sbrier", "get_nnem_brier_score", "estep",
    "marginal_density_yis", "density_weibull", "f_td", "sigma_est", "MSE_y_in",
    "CE" , "ibs_est_train", "ibs_est_test" ,"ISE_in" , "ISE_out" , "MSE_y_out" , 
    "MSE_y_in", "fixef_est", "scale_est", "shape_est", "D_est", "pij", "pij_true", "n"
  ),
  .combine = c, 
  .verbose = TRUE
) %dopar% {
  
  # Process each iteration in the batch sequentially
    # Your original computation here
    # Add debugging
  # iter <- 1
  c_rate <- 0.05
    cat("Processing iter:", iter, "\n")
     wkdir <- paste0("NNEM/Example2/k=2_Results_n=", n, "_c", c_rate)
file <- paste0(wkdir,"/nnem_size", iter, n,".RData")
load(file)
  theta <- result$theta
  # print("Error Variance:");theta$sigma_e;theta_orig$sigma_e;(bias_sigma <- abs(theta$sigma_e - theta_orig$sigma_e))
  if(!is.null(theta)){
    #  if(theta$sigma_e>0){
  sigma_est[iter] <- theta$sigma_e^2
  fixef_est[iter,] <- matrix(theta$betas,nrow=1, byrow=T)
  scale_est[iter,] <- matrix(theta$scale, nrow = 1, byrow = T) 
  shape_est[iter,] <- matrix(theta$shape, nrow=1, byrow = T)
  # omega_est[iter,] <- matrix(sqrt((theta$omegas)^2), nrow = 1, byrow=T)
  D_est[iter] <- theta$D
  pij[iter] <- list(theta$pj)
  
  pij_true[iter] <- list(result$theta_orig$pj)
  true_components<-numeric(n)
    
  for(j in seq_len(n)){
  red_c <- result$data_long[result$data_long$ID==j,]
  true_components[j] <- unique(red_c$g)
  }
# iter <- 1
# library(dplyr)
  # pi_true <- result$theta_orig$pj
      pred_y <- matrix(0, nrow(result$data_long), 2)
        Z <- model.matrix(~1, data = result$data_long)
    X <- model.matrix(~measure_time,data=result$data_long)

data <- result$data_long
data_wide <- result$data_wide
f_td <-  function(shape, scale, t, d) {
    ((dweibull(t,shape=shape, scale = scale)^d)) * (pweibull(t, shape = shape, scale =scale,lower.tail = F)^((1-d))) 
  }
comps <- estep(prior_mixing = theta$pj, theta)
comp.p <- comps$P
n_timepoints <- table(result$data_long$ID)     
bi <- list()
#? need to recompute class specific random effects 
# dgs <- tcrossprod(cholesky)
m <- length(theta$shape)
# omegas_tmp <- c(omegas,1)
q <- nrow(theta$D)
for(j in 1:m){
  #  j <- 1
    sigma_g <- result$theta$sigma_e # error variances 

    D_g <- theta$D
     # random effect variance-covariance matrices
    beta_g <- theta$betas[,j] # fixed effects

  #  iter<-1
  bigs<-matrix(nrow = n, ncol = q)
    for(i in 1:n){
    #   i<-2
      n_i <- n_timepoints[i]
      # k <- 1
      sim <- data
      subset<-sim[sim$ID==i,]
      Z_i <- Z[sim$ID==i, , drop = FALSE]
      X_i <- X[sim$ID==i, ,drop=FALSE]
       p<-ncol(X_i)
      Xi_beta <- X_i %*% beta_g
      # xbs[i]<-list(Xi_beta)
      # Z_i<-
      q<-ncol(Z_i)
      Y_i <- subset$y
      # ys<-list(Y_i)
      sigma<-sigma_g

      # 3by2 %*% 
      Vi <- as.matrix(Z_i%*% D_g%*% t(Z_i) + sigma_g^2 * diag(n_i))
      big <- D_g %*% t(Z_i) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
      bigs[i] <- big
      }
  bi[j] <- list(bigs)
}

    for(g in 1:2){

    mu <- NULL
    # g <- result$C[i]

    beta_tmp <-  theta$betas[,g]
    ranef_tmp <- bi[[g]]
  for (i in 1:n) {
    #  i <- 1
      # pi_tmp <- pij[[i]]
  # X_tmp <- cbind(1,result$data_long$measure_time)




    idx <- which(result$data_long$ID == i)
    Z_i <- Z[idx, , drop = FALSE]
    # V_i <- ((Z_i %*% (D) %*% t(Z_i)) + ((sd_e^2) * diag(length(idx))))
    X_i <- X[idx, ,drop=FALSE]
    XBg <- (X_i %*% beta_tmp + Z_i%*%ranef_tmp[i])
    # V_blocks[[i]] <- V_i
    mu <- c(mu, XBg)
  }
  pred_y[,g] <- mu
}

preds <- rowSums(pred_y * comp.p[result$data_long$ID,] )
# preds

  # for(k in 1:2){
  #   # k <- 1
  # pred_y[,k]<-  comp.p[result$data_long$ID,k] * ((X %*% theta$betas[,k]) + theta$ranef[[k]][result$data_long$ID])
  # }
  MSE_y_in[iter] <- Metrics::mse(actual = (result$data_long$y), predicted = (preds))
  
  CE[iter] <- mean(result$C  != true_components)
    data <- result$data_long
data_wide <- result$data_wide
result$theta_orig$betas <- t(result$theta_orig$betas)
true_comp <- estep(result$theta_orig$pj, result$theta_orig)$P
ISE_in[i]<- get_nnem_ISE(result, shapes = shape_est[i,], scales =scale_est[i,], mix.prob = comp.p, result$theta_orig, true_comp)

  ibs_est_train[iter]<- get_nnem_brier_score(shape_params=shape_est[iter,], scale_params=scale_est[iter,], mix.prob = comp.p, test_data=result$data_wide)$IBS$ibs
  

  # pij[i] <- list(theta$pj)
  # pij_true[i] <- list(pi_true)
  # est_surv <- pred_surv_marg(shape_params=shape_est[i,], scale_params=scale_est[i,], mix.prob = pij[[i]], data_wide=result$data_wide)

  #  times <- result$data_wide$event_time
  #  deltas <- result$data_wide$delta
  #  itx <- Surv(times, deltas)

  # ibs_est_train[i] <- get_nnem_brier_score(shape_params=shape_est[i,], scale_params=scale_est[i,], mix.prob = pij[[i]], test_data=result$data_wide)
   wkdir <- paste0("NNEM/Example2/k=2_Results_n=", n, "_c", c_rate)
    if(iter==1000){
      file <- paste0(wkdir,"/nnem_size", 1, n,".RData")}
    else{
      file <- paste0(wkdir,"/nnem_size", iter+1, n,".RData")
    }
  load(file)
  if(!is.null(result)){
data <- result$data_long
data_wide <- result$data_wide
result$theta_orig$betas <- t(result$theta_orig$betas)
true_comp <- estep(result$theta_orig$pj, result$theta_orig)$P
comp.p <- estep(result$theta_orig$pj, theta)$P

  ISE_out[i] <- get_nnem_ISE(result, shapes = length(shape_est[i,]), scales=scale_est[i,], mix.prob = comp.p,result$theta_orig, true_comp)
      pred_y <- matrix(0, nrow(result$data_long), 2)
        Z <- model.matrix(~1, data = result$data_long)
    X <- model.matrix(~measure_time,data=result$data_long)
    true_comp <- estep(result$theta$pj, theta)$P


    n_timepoints <- table(result$data_long$ID)     

  bi <- list()
#? need to recompute class specific random effects 
# dgs <- tcrossprod(cholesky)
m <- length(theta$shape)
# omegas_tmp <- c(omegas,1)
q <- nrow(theta$D)
for(j in 1:m){
  #  j <- 1
    sigma_g <- theta$sigma_e # error variances 

    D_g <- theta$D
     # random effect variance-covariance matrices
    beta_g <- theta$betas[,j] # fixed effects

  #  i<-1
  bigs<-matrix(nrow = n, ncol = q)
    for(i in 1:n){
    #   i<-2
      n_i <- n_timepoints[i]
      # k <- 1
      sim <- data
      subset<-sim[sim$ID==i,]
      Z_i <- Z[sim$ID==i, , drop = FALSE]
      X_i <- X[sim$ID==i, ,drop=FALSE]
       p<-ncol(X_i)
      Xi_beta <- X_i %*% beta_g
      # xbs[i]<-list(Xi_beta)
      # Z_i<-
      q<-ncol(Z_i)
      Y_i <- subset$y
      # ys<-list(Y_i)
      sigma<-sigma_g

      # 3by2 %*% 
      Vi <- as.matrix(Z_i%*% D_g%*% t(Z_i) + sigma_g^2 * diag(n_i))
      big <- D_g %*% t(Z_i) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
      bigs[i] <- big
      }
  bi[j] <- list(bigs)
}

    for(g in 1:2){

    mu <- NULL
    # g <- result$C[i]

    beta_tmp <-  theta$betas[,g]
    ranef_tmp <- bi[[g]]
  for (i in 1:n) {
    #  i <- 1
      # pi_tmp <- pij[[i]]
  # X_tmp <- cbind(1,result$data_long$measure_time)




    idx <- which(result$data_long$ID == i)
    Z_i <- Z[idx, , drop = FALSE]
    # V_i <- ((Z_i %*% (D) %*% t(Z_i)) + ((sd_e^2) * diag(length(idx))))
    X_i <- X[idx, ,drop=FALSE]
    XBg <- (X_i %*% beta_tmp + Z_i%*%ranef_tmp[i])
    # V_blocks[[i]] <- V_i
    mu <- c(mu, XBg)
  }
  pred_y[,g] <- mu
}

preds <- rowSums(pred_y * true_comp[result$data_long$ID,] )
preds


  MSE_y_out[iter] <- Metrics::mse(actual = (result$data_long$y), predicted = (preds))
  ibs_est_test[iter] <- get_nnem_brier_score(shape_params=shape_est[iter,], scale_params=scale_est[iter,], mix.prob = true_comp, test_data=result$data_wide)$IBS$ibs
  # est_surv <- pred_surv_marg(shape_params=shape_est[i,], scale_params=scale_est[i,], mix.prob = pij[[1]], data_wide=result$data_wide)
  # times <- result$data_wide$event_time
  # deltas <- result$data_wide$delta
# itx <- Surv(times, deltas)
  # est_ibs <- sbrier(obj=itx,pred=est_surv)
      # ibs_est_test[i] <- sbrier(obj = itx, est_surv)


  }
  }
  # pred_c[i] <- list(result$C)
  
  # extra <- nchar('||100%')
  # width <- options()$width
  # iter <- iter
  # maxit <- B
  # step <- round(iter / maxit * (width - extra))
  # text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
  #                 strrep(' ', width - step - extra), round(iter / maxit * 100))
  # cat("\nSimulations Data Offload Progress:\n",text, '\n')
  # # Sys.sleep(0.5)
  # cat(if (iter == maxit) '\n' else '\014')
        list(list(
          iter=iter,
      sigma_est =  sigma_est[iter],
      fixef_est= fixef_est[iter,], 
  scale_est = scale_est[iter,],
  shape_est = shape_est[iter,],
  # omega_est[iter,],
  D_est=D_est[iter],
  pij=pij[[iter]],
  MSE_y_in= MSE_y_in[iter],
  MSE_y_out= MSE_y_out[iter],
  ibs_est_train=ibs_est_train[iter],
  ibs_est_test =ibs_est_test[iter],
   CE=CE[iter]
      # oob_idx = oob_idx,
      # oob_pred = oob_class
    ))
}
    for(result in offload) {
  if(!is.null(result)) {
    (iter <- result$iter)
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
}
    results_filename <- paste0("results_n", n, "_c", gsub("\\.", "", c_rate), ".RData")
    save(sigma_est, fixef_est, scale_est, shape_est, D_est, pij, pij_true, 
         CE, ibs_est_train, ibs_est_test, ISE_in, ISE_out, MSE_y_out, MSE_y_in,
         file = results_filename)
    
    cat("Completed n =", n, ", censoring rate =", c_rate, "\n")
  }
}
toc()




  class_err <- mean(CE)


# est_c <- numeric(200)
# for(i in seq_len(B)){
# est_c[i] <- mean(c_true[[i]]!=pred_c[[i]])
# }
# est_surv <- pred_surv_marg(shape_params=shape_est[1,], scale_params=scale_est[1,], mix.prob = pij[[1]], data_wide=result$data_wide)
# times <- result$data_wide$event_time
# deltas <- result$data_wide$delta
# itx <- Surv(times, deltas)
# est_ibs <- sbrier(obj=itx,pred=est_surv)
class_err <- mean(CE)
# idxs_to_keep <- which(ibs_est_test>0) #& scale_est[,2]<(50+3*sd(scale_est[,2])) & scale_est[,2]>50-3*sd(scale_est[,2]))
 #& scale_est[,2]<(50+3*sd(scale_est[,2])) & scale_est[,2]>50-3*sd(scale_est[,2]))
# ibs_est_test_est<- mean(ibs_est_test[idxs_to_keep])
idxs_to_keep <- which(ibs_est_train>0)
ibs_est_train_est <- mean(ibs_est_train[idxs_to_keep])

ise_out <- mean(ISE_out[ISE_out>0])
ise_in <- mean(ISE_in[ISE_in>0])
theta_orig <- result$theta_orig
theta_orig$betas <- c(t(theta_orig$betas))
idxs_to_keep <- which(sigma_est>0) #& scale_est[,2]<(50+3*sd(scale_est[,2])) & scale_est[,2]>50-3*sd(scale_est[,2]))
sigma_est <- sigma_est[idxs_to_keep]
D_est <- D_est[idxs_to_keep]
fixef_est <- fixef_est[idxs_to_keep,]
scale_est <- scale_est[idxs_to_keep,]
shape_est <- shape_est[idxs_to_keep,]
cat("\nError Variance:",mean(sigma_est),theta_orig$sigma_e^2)#(bias_sigma <- abs(theta$sigma_e - theta_orig$sigma_e));mean(sigma_est)
# print("Cholesky");theta_orig$cholesky; colMeans(chol_est); #(bias_cholesky <- abs(theta$cholesky - theta_orig$cholesky));colMeans(chol_est)
cat("\nRE Variance-Covariance",mean(D_est), theta_orig$D[1]);#(bias_D <- abs(vech(theta$D)-vech(theta_orig$D)));colMeans(D_est)
# print("Omegas:");colMeans(sqrt(omega_est^2)); theta_orig$omegas;#(bias_omegas <- abs(sqrt((theta$omegas)^2) - theta_orig$omegas));colMeans(omega_est)
cat("\nLongitudinal Fixef",colMeans(fixef_est), theta_orig$betas);#(bias_beta <- abs(theta$betas - theta_orig$betas));colMeans(fixef_est)
cat("\nWeibull Scales",colMeans(scale_est), theta_orig$scale);#(bias_scale <- abs((theta$scale)-(theta_orig$scale)));colMeans(scale_est)
cat("\nWeibull Shapes", colMeans(shape_est), theta_orig$shape)# (bias_shape<- abs(theta$shape - theta_orig$shape));colMeans(shape_est)

GOF <- rbind(class_err,ibs_est_test_est,ibs_est_train_est, ise_out, ise_in)

bias_sigma <- abs(mean(sqrt(sigma_est))-theta_orig$sigma_e);bias_sigma;
bias_scale <- abs(colMeans(scale_est)-matrix(theta_orig$scale, nrow = 1, byrow = T));bias_scale;colMeans(scale_est)
bias_shape <- abs(colMeans(shape_est)-matrix(theta_orig$shape, nrow = 1, byrow = T));bias_shape;colMeans(shape_est)
bias_fixef <- abs(colMeans(fixef_est)-matrix(theta_orig$betas, nrow = 1, byrow = T));bias_fixef; MAB_fixef <- (1/(3*(2+1)))*sum(bias_fixef); MAB_fixef;colMeans(fixef_est)
bias_D <- abs(mean(D_est)-matrix(theta_orig$D, nrow = 1, byrow = T));bias_D;mean(D_est)
# bias_cholesky <- abs(colMeans(chol_est)-matrix(theta_orig$cholesky[lower.tri(theta$cholesky, diag = T)], nrow = 1, byrow = T));bias_cholesky;colMeans(chol_est)
# bias_omegas <- abs(colMeans(omega_est)-matrix(theta_orig$omegas, nrow = 1, byrow = T));bias_omegas;colMeans(omega_est)


# install.packages("caret")
# MSE
B <- length(sigma_est)
true_p <- MSE_p <- numeric(2)
for(k in seq_len(2)){
MAB_pi <- true_pi <- matrix(ncol=B, nrow = n)
for(i in seq_len(B)){
  if(!is.null(pij[i][[1]])){
  # print(k)
MAB_pi[,i] <- pij[i][[1]][,k]
true_pi[,i] <- pij_true[i][[1]][,k]
}}
pi2 <- mean(abs(rowMeans(MAB_pi, na.rm = T) - rowMeans(true_pi, na.rm = T)));
pi3 <- mean(rowMeans(((MAB_pi) - (true_pi))^2, na.rm = T))
true_p[k] <- c(pi2)
MSE_p[k] <- pi3

}
mix.prop.bias <- true_p

MSE_p
true_p
matE <- rep(theta_orig$sigma_e,B)
MSE_sigma <- 1/B * t(sqrt(sigma_est)-matE)%*%(sqrt(sigma_est)-matE)


matB <- matrix(rep(theta_orig$betas,B), nrow=B, ncol=4, byrow = T)
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


matB <- matrix(rep(theta_orig$betas,B), nrow=B, ncol=4, byrow = T)
# correction factor
MSE_beta <- sapply(1:4, function(i) Metrics::mse(matB[,i], fixef_est[,i]))
MAB_beta <-   colMeans(abs(fixef_est-matB)); MAB_beta; 
MAB_beta <- sapply(1:4, function(i) Metrics::mdae(matB[,i], fixef_est[,i]))

#Tau (Cholesky of D)
# matT <- matrix(rep(theta_orig$cholesky[lower.tri(theta_orig$cholesky, diag=T)],B), nrow=B, ncol=3, byrow = T)
# ave_T <- colMeans(chol_est)
# bias_cholesky <- abs(ave_T-matrix(theta_orig$cholesky[lower.tri(theta$cholesky, diag = T)], nrow = 1, byrow = T));bias_cholesky;
# MSE_T <- colMeans((chol_est-matT)^2); sqrt(MSE_T)
# MAB_T <-   colMeans(abs(chol_est-matT)); MAB_T

#D
bias_D <- abs(mean(sqrt(D_est))-matrix(sqrt(theta_orig$D[1]), nrow = 1, byrow = T));bias_D;
ave_D <- mean(sqrt(D_est))
matD <- matrix(theta_orig$D, nrow = B, ncol=1, byrow = T)
# MSE_D <- mean((sqrt(D_est)-sqrt(matD))^2); (MSE_D)
# (MSE_D <- var(sqrt(D_est))+(bias_D)^2)
MSE_D <- sapply(1:1, function(i) Metrics::mse(sqrt(matD), sqrt(D_est)))
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
    mean(sqrt(sigma_est)),
    colMeans(fixef_est),
    mean(sqrt(D_est)),
    colMeans(shape_est),
    colMeans(scale_est)
    # colMeans(omega_est)
  ),
  
  Bias = c(
    # Bias calculations
    true_p,
    bias_sigma,
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
                       caption = "Bias and MSE Results for NNEM K=2 Parameter Estimates (B=1000 simulations, 5% censoring)",
                       label = "tab:bias_mse_results",
                       digits = 3)

# Print LaTeX code
print(xtable_output, 
      type = "latex",
      include.rownames = FALSE,
      sanitize.text.function = function(x) x,  # Preserve LaTeX math notation
      caption.placement = "top",
      table.placement = "H")

print("=== XTABLE OUTPUT ===")
xtable_output <- xtable(GOF, 
                       caption = paste0("Metrics of NNEM for sample size ",n),
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


