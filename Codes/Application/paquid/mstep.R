#! M - step
gform=~CEP+male
mstep <- function(postprob, theta, data_wide, data_long, ranform, long_form_g, long_form_o, gform, survform){
# theta <- list(scale=scale_params, shape=shape_params, gammas=gammas_matrix, betas=betas_matrix, D=D_matrices, sigma_e=sigma)

# postprob <- estep.result$P
  # postprob <- e_step_result$P
  n <- nrow(data_wide)
  m <- nrow(betas_matrix)
  x <- model.frame(gform, data_wide)
# mixing proportions using nnet

    # For more than 2 components, use multinomial
    # postprob <- e_step_result$postprob

    nn_model <- nnet(
      x = x,          # Input features
      y = postprob,          # Target is component probabilities
      size = 10,      # Hidden layer size
      decay = 0.01,   #! Weight decay: this is something which can be tweaked
      maxit = 5000,   # Max iterations
      trace = FALSE,  # No output during training
      softmax = TRUE,  # Use softmax for multinomial
      entropy=TRUE # maximum conditional likelihood fitting
    )
    
    # Predict probabilities (this code doesn't work for some stupid reason, 
  # extract directly from the nnet object)
    # prob <- predict(nn_model, data_wide$X3)
    new_prior <- fitted(nn_model)


print("New priors")
print(head(new_prior))
  
random_effects <- theta$ranef

#? compute random effects matrices (D)
# print(theta$cholesky)
# init_params_j <- matrix(c(sqrt(theta$omegas),theta$cholesky[lower.tri(theta$cholesky, diag = T)]))

# neg_loglik_omegas_Dg <- function(parms_j, data, ranef, prior, posterior, m=3){

#   cholesky <- metafor::vec2mat(parms_j[(m):(m+2)], diag=T)
#   cholesky[1,2] <- 0
#   # D_g <- diag(2)
#   print(cholesky)

#   omegas <- parms_j[1:(m-1)]

#   print(omegas)


# random_effects_density <- function(u, G, omega, cholesky_parameterization = FALSE, log = FALSE) {
  
#   if (cholesky_parameterization) {
#     # G is provided as Cholesky factor L, so G_actual = L %*% t(L)
#     L <- G
#     G_actual <- (omega^2)*(L %*% t(L))
#   } else {
#     # G is the actual covariance matrix
#     L <- G
#     G_actual <- (omega^2)*(L %*% t(L))
#     # G_actual <- G
#   }
  
#   # f(u) = (2π)^(-q/2) |G|^(-1/2) exp(-1/2 u' G^(-1) u)
  
#   if (log) {
#     q <- length(u)
    
#     if (cholesky_parameterization) {
#       # Efficient computation using Cholesky factor
#       # log|G| = 2 * sum(log(diag(L)))
#       log_det_G <- 2 * sum(log(diag(t(chol(G_actual)))))
      
#       # Solve L w = u for w, then u' G^(-1) u = w' w
#       w <- solve(t(chol(G_actual)), u)
#       quadratic_form <- t(w)%*%w
      
#     } else {
#       # Standard computation
#       log_det_G <- determinant(G_actual, logarithm = TRUE)$modulus[1]
#       quadratic_form <- t(u) %*% solve(G_actual) %*% u
#     }
    
#     log_constant <- -q/2 * log(2 * pi) - 0.5 * log_det_G
#     log_quadratic <- -0.5 * quadratic_form
    
#     return(log_constant + log_quadratic)
    
#   } else {
#     # Use dmvnorm for numerical stability
#     return(dmvnorm(u, mean = rep(0, length(u)), sigma = G_actual, log = T))
#   }
# }
# # ll <- matrix(0,nrow=100, ncol=m)
# ll <- 0
# omegas <- c(omegas,1)
# for(i in 1:nrow(ranef[[1]])){
#   for(k in 1:m){
#     ll <- ll + posterior[i,k]*(log(prior[i,k]+0.00001)+random_effects_density(u=ranef[[k]][i,],G=cholesky, cholesky_parameterization = T, log = T, omega = omegas[k]))
#   }
# }
#   # ll <- sum(rowSums(ll))
#   print(-ll)
#   return(-ll)
# }
# # this is increases the log-likelihood (nll should be positive)
# neg_loglik_omegas_Dg(init_params_j,data=data,  ranef=random_effects, prior=new_prior, posterior=postprob)


# library('optimParallel')
# n_cores <- detectCores() - 1
# cl <- makeCluster(n_cores)     # set the number of processor cores
# setDefaultCluster(cl=cl)
# clusterExport(cl,"dmvnorm", envir = environment())
# library(tictoc)
# tic()
# opt_result <- optimParallel(as.vector(init_params_j), neg_loglik_omegas_Dg, ranef=random_effects, prior=new_prior, posterior=postprob,hessian = T, parallel = list(cl=cl, loginfo=T), data=data, verbose=T)
# toc()
# stopCluster(cl)

# cholesky <- metafor::vec2mat(opt_result$par[(m):(m+2)], diag=T)/(1)
# cholesky[1,2] <- 0
# print(tcrossprod(cholesky))
# # D_g <- matrix(c(parms_j[4], parms_j[5], parms_j[5], parms_j[6]), nrow = T)
# omegas <- opt_result$par[1:(m-1)]
# # print(cholesky)
#     if(!is.null(opt_result$hessian)& !is.singular.matrix(opt_result$hessian)){
# sds <- sqrt(diag(matrix.inverse(opt_result$hessian)))

# chol_sd <- metafor::vec2mat(sds[m:(m+2)], diag=T)
# chol_sd[1,2] <- 0
# omegas_sd <- sds[1:(m-1)]
#   } else{
#       omegas_sd <- NULL
#       chol_sd <- NULL
#     }


# # #   dgs_sds[[k]] <- metafor::vec2mat(sqrt(diag(matrix.inverse(opt_result$hessian))), diag = T)

# #   # 
# # # neg_loglik_omegas_Dg(init_params_j,data=data,  ranef=random_effects, prior=new_prior, posterior=postprob)
# # cat("Contrib. Omegas:",-neg_loglik_omegas_Dg(matrix(c(sqrt(omegas),cholesky[lower.tri(theta$cholesky, diag = T)])),data=data,  ranef=random_effects, prior=new_prior, posterior=postprob))




# log_sde <- log(theta$sigma_e^2); log_sde

# #? compute group specific betas and error variance (this takes an age to compute)
# #' Alternative: Using individual normal densities (more stable)

# #' Total conditional negative log-likelihood (sum over all individuals)
# total_nll_betas <- function(parms, y, random_effects, X, Z, postprob, prior, data) {

# lli_stable <- function(y_i, u_i, X_i, Z_i, beta, sigma_e) {
  
#   # Conditional mean
#   mu_i <- X_i %*% beta + Z_i %*% u_i
  
#   ll_i <- mvnfast::dmvn(t(y_i), mu_i,(sigma_e)*diag(length(mu_i)), log = T)
#   # Sum of negative log-densities
#   # ll_i <- sum(dnorm(y_i, mean = mu_i, sd = sigma_e, log = T))
  
#   return(ll_i)
# }
#     betas_matrix <- matrix(parms[-7],ncol = 2) # first p*m parameters are the betas
#   print(betas_matrix)
#   sd_e <- exp(parms[7]) #the p*m+1 parameter should be the standard error
#   print(sd_e)
#   individual_indicators <- data$ID
#   n_subjects <- max(individual_indicators)
#   m <- nrow(betas_matrix)

#   loglik <- matrix(0, n_subjects, m)
#   ll <- 0
  

  
#   for (i in 1:n_subjects) {
#       for(k in 1:m){
#   u_matrix <- random_effects[[k]]
#     # Extract data for individual i
#     idx_i <- which(individual_indicators == i)
#     y_i <- y[idx_i]
#     X_i <- X[idx_i,]
#     Z_i <- Z[idx_i,]
#     u_i <- u_matrix[i, ]  # ith row = [u₀ᵢ, u₁ᵢ]
    
#     # Compute likelihood for individual i
#     ll <- ll + (postprob[i,k]+0.00001)*(lli_stable(y_i, u_i = u_i, X_i, Z_i, beta =betas_matrix[k,], sigma_e = sd_e))

#     # nll_i <- lli_stable(y_i, u_i, X_i, Z_i, beta, sigma_e)
#     # total_nll <- total_nll + nll_i
#   }}
#   # total_nll <- sum(rowSums(loglik))
  
#   return(-ll)
# }
# init_params_j <- c(theta$betas, log_sde)
# cat("\nPrevious Contribution for Betas\n",-total_nll_betas(init_params_j,data=data,postprob=postprob,prior=pi_init, random_effects=random_effects, y=y, X=X, Z=Z))
# n_cores <- detectCores() - 1
# cl <- makeCluster(n_cores)     # set the number of processor cores
# setDefaultCluster(cl=cl)
# # opt_result <- optim(init_params_j, total_nll_betas, data=data, prior=pi_init, postprob=postprob, y=y, X=X, Z=Z, random_effects=random_effects, hessian = T, control = list(maxit=1000))
# opt_result <- optimParallel(init_params_j, total_nll_betas, prior=pi_init, postprob=postprob, y=y, X=X, Z=Z, random_effects=random_effects,data=data,hessian=T,parallel = list(cl=cl, loginfo=T))
# #   ?optim
# library(marqLevAlg)
  # opt_result <- mla(b=init_params_j, fn=total_nll_betas, prior=pi_init, postprob=postprob, y=y, X=X, Z=Z, random_effects=random_effects,data=data,nproc = n_cores)
  # hessian <- summary(opt_result, loglik=T)
  # init_params_j <- c(theta$betas, theta$sigma_e)
  # neg_loglik_betak_sde <- function(parms_j, data){
  #   # parms_j <- init_params_j
  #   betas_matrix <- matrix(parms_j[-7],ncol = 2) # first p*m parameters are the betas
  #   print(betas_matrix)
  #   sd_e <- parms_j[7] #the p*m+1 parameter should be the standard error
  #   print(sd_e)
  #   f_yCondb <- log_conditional_density_yis(data_long=data, betas_matrix = betas_matrix, sigmas=sd_e, random_effects = random_effects)
  #   mixed <- sum(rowSums(postprob*(f_yCondb)))
  #   return(-mixed)
  # }
  
#  opt_result <- optim(init_params_j, neg_loglik_betak_sde,hessian=T, data=data, method="BFGS")
  
  #  loglik_betak_sde <- function(parms_j, data){
  # #   # parms_j <- init_params_j
  # #   betas_matrix <- matrix(parms_j[-9],ncol = 2) # first p*m parameters are the betas
  # #   # print(betas_matrix)
  # #   sd_e <- parms_j[9] #the p*m+1 parameter should be the standard error
  # #   # print(sd_e)
  # #   f_yCondb <- log_conditional_density_yis(data_long=data_long, betas_matrix = betas_matrix, sigmas=sd_e, random_effects = ranef)
  # #   ll_betak_sde <- sum(rowSums(posteriors*(f_yCondb)))
  # #   return(mixed)
  # # }
  
  # cl <- makeCluster(n_cores)     # set the number of processor cores
  # setDefaultCluster(cl=cl)
  # clusterExport(cl, c(#bis
  #   "log_conditional_density_yis", "random_effects", "nis", "postprob",
  #   # betas_sds
  #   "m", "p", "q"), envir = environment())
  # tic()
  # opt_result <- tryCatch({optimParallel(as.vector(init_params_j), neg_loglik_betak_sde, hessian = T, data=data, parallel = list(cl=cl, loginfo=T), verbose=T)},error=function(e){{optimParallel(as.vector(init_params_j), neg_loglik_betak_sde, data=data, parallel = list(cl=cl, loginfo=T), verbose=T)}})
  # toc()
  # stopCluster(cl)

init <- c(theta$betas, theta$beta_o, theta$cholesky[lower.tri(theta$cholesky, diag = T)],log(theta$sigma_e)) 
# parms <- init
# data <- data_long
# posterior <- postprob
# theta_true <- theta
nll_marg <- function(parms, data, posterior, long_form_g, long_form_o, ranform, theta_true){
  
  l1 <<- nrow(theta_true$betas)*ncol(theta_true$betas)
l2 <<- sum(lower.tri(theta_true$cholesky, diag = T))
log_marginal_density_yis <- function(data_long, theta, long_form_g, long_form_o,ranform) {
  # parms prep
  betas_matrix <- theta$betas
  sigmas <- theta$sigma_e
  # D_matrices <- theta$D
  beta_o <- theta$beta_o # overall sample fixed effects parameters

  # data prep
  y <- model.frame(long_form_g,data_long)[,1]
  Xg <- model.matrix(long_form_g, data_long)
  Xo <- model.matrix(long_form_o, data_long)
  Zg <- model.matrix(ranform, data_long)
  data_long$y <- y

  N <- length(unique(data_long$ID));#print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix);#print(m);#print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  nis<-table(data_long$ID)
  IDs<-unique(data_long$ID)
  # length(nis)

  D_matrices <- tcrossprod(theta$cholesky);# print(D_matrices)
  for (k in 1:m) {
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects
  
    for(i in seq_len(length(IDs))){

      # i<-2
      # n_i <- nis[i]
      # k <- 1
      ids <- IDs[i] 
      
      # data_long <- data;i <- 1
      subset <- data_long[data_long$ID == ids, ] #! ensure that IDs are consistently named across datasets
            # subset<-sim[sim$ID==i,]
      Z_ig <- Zg[data_long$ID == ids, ,drop = FALSE]    # Random Effects Design Matrix
      X_ig <- Xg[data_long$ID==ids,,drop=FALSE ]# Cluster-specific Fixed Effects Design Matrix
      X_io <- Xo[data_long$ID==ids,,drop=FALSE] # Overall Fixed Effects Design Matrix (*no intercept)

      p <- ncol(X_ig)
      Xi_beta_k <- matrix(X_ig %*% beta_g + X_io %*% beta_o) # Fixed Effects Linear Predictor
   
      q <- ncol(Z_ig)
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      V_ig <- as.matrix(Z_ig %*% (D_g) %*% t(Z_ig) + sigma_g^2 * diag(nrow(X_ig)))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig, log = T)
    }
    }
  return(f_yi)
}


  init <- parms
  betas_matrix <- matrix(init[1:l1], nrow = nrow(theta_true$betas), ncol=ncol(theta_true$betas))
  beta_o <- init[l1+1]
  cholesky <- metafor::vec2mat(init[(l1+2):(l1+l2+1)], diag = T)
  # cholesky[3,1] <- cholesky[2,2]
  # cholesky[2,2] <- cholesky[1,3]

  cholesky[1,2] <- cholesky[1,3] <- cholesky[2,3]<- 0
  #print(tcrossprod(cholesky))
  sd_e <- exp(init[length(init)])
  #print(sd_e)
    tmp_theta <- list(betas=betas_matrix,beta_o=beta_o, sigma_e=sd_e,cholesky=cholesky)
ysu <- log_marginal_density_yis(
  data_long = data, 
theta=tmp_theta,
long_form_g = long_form_g,
long_form_o = long_form_o,
ranform = ranform)
  
  llk <- sum(rowSums(posterior*(ysu)))
  return(-llk)
}
# nll_marg(init, data, postprob, long_form_g, long_form_o, ranform, theta)
library(parallel)
library(marqLevAlg)
n_cores <- detectCores() - 1
opt_result <- marqLevAlg::mla(init, fn=nll_marg,data=datalong, posterior=postprob,nproc=n_cores, print.info=T, long_form_g=long_form_g, long_form_o=long_form_o, ranform=ranform, theta_true=theta)
# opt_result <- obj
#  cl <- makeCluster(n_cores)     # set the number of processor cores
#  setDefaultCluster(cl=cl)
# opt_result <- optimParallel(init, fn=nll_marg, data=data, posterior=postprob,hessian=T, parallel = list(cl=cl))
cat("============================================")
cat("\nNN-JLCM Longitudinal Parameter Estimates:\n")
cat("============================================")
  sd_info_fixef <- summary(opt_result, loglik=T,, digits = 10)

  l1 <- nrow(theta$betas)*ncol(theta$betas)
l2 <- sum(lower.tri(theta$cholesky, diag = T))
# obj <- opt_result
(betas_matrix <- matrix(opt_result$b[1:l1], nrow = m, ncol=3))
  beta_o <- opt_result$b[(l1+1)]
  cholesky <- metafor::vec2mat(opt_result$b[(l1+2):(l1+1+l2)], diag=T)
  cholesky[1,2] <- cholesky[1,3] <- cholesky[2,3]<- 0
   (D <- tcrossprod(cholesky))
  (sd_e <- exp(opt_result$b[length(opt_result$b)]))

# theta <- theta_orig
shapes <- as.matrix(theta$shape)
scales <- as.matrix(theta$scale)
# gammas_matrix <- theta$gammas
#? compute random effects matrices (D)

# cat("Variance Component Matrices and Associated Standard Errors:\n")
# print(chol_sd)
 
  # as.numeric(sd_info$coef)
# print(tcrossprod(cholesky),quote=F)


# print(omegas_sd)

# cat("Residual Standard Deviation (s.e):\n",sd_e," (",sd_info$SE.coef[19],") ", sep="")
# betas_matrix <- matrix(opt_result$par[-7],ncol=2)
# # # opt_result$par
# sd_e <- sqrt(exp(opt_result$par[7]))
# abs(init_pault$par[9]; 
# ?optim

# if()
# if(!is.null(opt_result$hessian)& !is.singular.matrix(opt_result$hessian)){
# betas_matrix_sds <- matrix(sqrt(diag(matrix.inverse(opt_result$hessian)))[-7], ncol=2); cat("\nBetas:",betas_matrix,"\nStandard Errors:",betas_matrix_sds)
# sd_e_sd <- sqrt(diag(matrix.inverse(opt_result$hessian)))[7]; cat("\nError Standard Deviation:", sd_e,"\nStandard Errors:",sd_e_sd)
# #! update names
# } else{
#   betas_matrix_sds <- NULL
#   sd_e_sd <- NULL
# }
# if(!is.null(opt_result$hessian)& !is.singular.matrix(opt_result$hessian)){
# betas_matrix_sds <- matrix(sqrt(diag(matrix.inverse(opt_result$hessian)))[-7], ncol=2); cat("\nBetas:",betas_matrix,"\nStandard Errors:",betas_matrix_sds)
# sd_e_sd <- sqrt(diag(matrix.inverse(opt_result$hessian)))[7]; cat("\nError Standard Deviation:", sd_e,"\nStandard Errors:",sd_e_sd)
# #! update names
# } else{
#   betas_matrix_sds <- NULL
#   sd_e_sd <- NULL
# }
# init_params_j <- c(betas_matrix, sd_e)

# cat("\nLog-likelihood Contribution - Betas and Sigma:\n",-total_nll_betas(init_params_j,data=data,postprob=postprob,prior=pi_init, random_effects=random_effects, y=y, X=X, Z=Z))
# # beta_matrix <- beta_gg
D_matrices <- D
# cholesky <- cholesky
sigmas <- sd_e
pi_init <- new_prior

# shapes <- as.matrix(theta$shape)
# scales <- as.matrix(theta$scale)
# gammas_matrix <- theta$gammas
# wx <- cbind(data_wide$X1, data_wide$X2)
# posterior <- c
# shapes_sd <- shapes
# scales_sd <- scales
# gammas_sd <- gammas_matrix
# wb <- density_weibull(shapes, scales)
# shapes
# scales
# gammas_matrix
parms <- c(log(shapes), log(scales), theta$gamma_matrix)
neg_loglik_wb <- function(parms_j, data_wide, posterior, ng){
  # parms_j <- parms
  log_shapes <- matrix(parms_j[1:m],ncol = 1); log_shapes # first m parameters are the shapes
  log_scales <- matrix(parms_j[(m+1):(2*m)]); log_scales # next m parameters are the scales
  gammas_matrix <- as.matrix(parms_j[(2*m+1):(length(parms_j))])

  f_td <-  function(shape, scale, t, d) {
    (d*(dweibull(t,shape=shape, scale = scale,log = T))) + (1-d)*(pweibull(t, shape = shape, scale =scale,lower.tail = F,log.p = T)) 
  }

  density_weibull <- function(shape_params, scale_params, gammas_params, data_wide, survform) {

    survform <- update(survform, ~.-1)
    W <- model.matrix(survform, data_wide)
    N <- nrow(data_wide)
    k <- nrow(as.matrix(shape_params))


    f_tds <- matrix(0,nrow = N, ncol = k)
    # w <- model.matrix(~ -1+X3+X4+X5, data=data_wide)
    delta <- data_wide$delta
    # initial density function for Weibull
    
    for (j in 1:k) {
      # j <- 1
      shape <- shape_params[j]
      scale <- scale_params[j]
      gammas <- gammas_params
      linear_pred <- exp((W%*%gammas)/shape) #updated from Prof. Breheny's lecture slides p18
      # linear_pred <- exp()
      delta <- delta
      new_scales <- (linear_pred)*scale
      # fixed effects in survival component
      # wxi <- as.matrix(w)
      ti <- data_wide$event_time
      for(i in 1:N){
      f_tds[i, j] <- f_td(shape, new_scales[i], ti[i], delta[i])
      }
      }
    return(f_tds)
  }
  scales <- exp(log_scales)
  shapes <- exp(log_shapes)
  lwb <-  density_weibull(shape_params = shapes, scale_params = scales, gammas_params = gammas_matrix, data_wide=data_wide,survform)
  component_likelihood <- sum(rowSums(posterior * (lwb)))
  return(-component_likelihood)
}
neg_loglik_wb(parms, data_wide = data_wide, posterior = postprob, ng=m)
opt_result_wb <- tryCatch({
  marqLevAlg::mla(b=parms, fn=neg_loglik_wb, data_wide=data_wide, posterior=postprob, ng=4, nproc=n_cores, print.info=T)}, error = function(e) {
    optim(par=parms, fn=neg_loglik_wb, data_wide=data_wide, posterior=postprob)})

cat("==========================================")
cat("\nNN-JLCM Survival Parameter Estimates:\n")
cat("==========================================")
sd_res_Td <- summary(opt_result_wb, loglik=T)
# sd_res_Td[,1:3] <- as.numeric(unlist(sd_res_Td[,1:3]))
# sd_res_Td[1:(2*m),c(1,5:6)] <- exp(sd_res_Td[1:(2*m),c(1,5:6)])
shapes <- exp(matrix(opt_result_wb$b[1:m],ncol = 1)) # first m bameters are the shapes
# sd_res_Td
scales <- exp(matrix(opt_result_wb$b[(m+1):(2*m)])) # next m bameters are the scales
gammas_matrix <- as.matrix(opt_result_wb$b[(2*m+1):(length(parms))])
# print("Successful convergeance of Weibull submodel parameters")
# std_errs <- sqrt(diag(matrix.inverse(opt_result$)))
#   if(!is.null(opt_result$hessian)& !is.singular.matrix(opt_result$hessian)){
#     std_errs <- sqrt(diag(matrix.inverse(opt_result$hessian)))
# # std_errs <- as.numeric(summary(opt_result, loglik=TRUE)$SE.coef)
# shapes_sd <- std_errs[1:m]
# scales_sd <- std_errs[(m+1):(2*m)]
#   }else{
   shapes_sd <- NULL
   scales_sd <- NULL
#   }

# cat("\nWeibull Shapes:", shapes, "\nStandard Errors:", shapes_sd)
# cat("\nWeibull Scales:", scales, "\nStandard Errors:", scales_sd)

# parms <- c(log(shapes), log(scales))
# cat("\nContribution for Weibull:", -neg_loglik_wb(parms, data_wide = data_wide, posterior = postprob))
# -neg_loglik_wb(parms, data_wide = data_wide, posterior = postprob)
  bi <- list()
  Xg <- model.matrix(long_form_g, data=datalong)
Xo <- model.matrix(long_form_o, datalong)
Zg <-  model.matrix(ranform, data=datalong)
for(j in 1:m){
# j <- 1
    sigma_g <- sd_e # error variances 

    D_g <- D
     # random effect variance-covariance matrices
    beta_g <- betas_matrix[j,] # fixed effects
    beta_o <- beta_o; print(beta_o)

  # i<-1
  bigs<-matrix(nrow = n, ncol = q)
  seq <- data_wide$ID

    for(i in seq_len(length(seq))){

      # i<-2
      # n_i <- nis[i]
      # k <- 1
      ids <- seq[i]
      sim <- datalong
      subset<-sim[sim$ID==ids,]
      Z_ig <- Zg[sim$ID==ids, , drop = FALSE]
      X_ig <- Xg[sim$ID==ids, ,drop=FALSE]
      X_io <- Xo[sim$ID==ids, ,drop=FALSE]

      p<-ncol(X_ig)
      Xi_beta <- X_ig%*%beta_g + X_io%*%beta_o
      # xbs[i]<-list(Xi_beta)
      q<-ncol(Z_ig)
      Y_i <- subset$y
      # ys<-list(Y_i)
      sigma<-sigma_g
      n_i <- length(subset$ID)
      # 3by2 %*% 
      Vi <- as.matrix(Z_ig%*% D_g%*% t(Z_ig) + sigma_g^2 * diag(n_i))
      big <- D_g %*% t(Z_ig) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
      bigs[i,] <- big
      }
  bi[j] <- list(bigs)
}


# standard_errors <- list()
# update parameters
 theta <- list(pj=new_prior, scale=scales, shape=shapes, betas=betas_matrix, beta_o=beta_o, gamma_matrix=gammas_matrix, ranef=bi, D=D, cholesky=cholesky, sigma_e=sd_e)
  results_standard_errors <- list("Longitudinal Submodel"=sd_info_fixef, "Survival Submodel"=sd_res_Td)
 ret <- list("Parameters"=theta, "Model Estimates & S.E."=results_standard_errors)
return(ret)
}


# pj <- estep.result$P
# # # library(tictoc)
# # # tic()
# # library(matrixcalc)
# m_step_result <- mstep(postprob=pj,theta, data_wide, data_long, ranform, long_form_g, long_form_o, gform, survform)
# cat("\n")
# toc()