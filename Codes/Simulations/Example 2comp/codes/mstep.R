#! M - step

mstep <- function(postprob, theta){
# theta <- list(scale=scale_params, shape=shape_params, gammas=gammas_matrix, betas=betas_matrix, D=D_matrices, sigma_e=sigma)

   postprob <- estep.result$P
  
  m <- nrow(betas_matrix)
# mixing proportions using nnet
  if (m == 2) {
    # Train on first component probability only
nn_model <- nnet(
    x = x, # Input features
      y = postprob,      # Target is first component probability
      size = 10,       # Hidden layer size
      decay = 0.01,    # Weight decay
      maxit = 5000,    # Max iterations
      trace = FALSE,  # No output during training
      softmax = TRUE,  # Use softmax for multinomial
      entropy=TRUE # maximum conditional likelihood fitting
    )
    
    # Predict probabilities
    prob <- fitted(nn_model)
    
    # Ensure values are between 0 and 1
    # prob1 <- pmin(pmax(prob1, 0), 1)
    
    # # Create complete probability matrix
    # prob <- cbind(prob1, 1 - prob1)
    new_prior <- prob
  } else {
    # For more than 2 components, use multinomial
    # postprob <- e_step_result$postprob

    nn_model <- nnet(
      x = x,          # Input features
      y = postprob,          # Target is component probabilities
      size = 5,      # Hidden layer size
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
}

print("New priors")
print(head(new_prior))
  
random_effects <- theta$ranef


#? compute random effects matrices (D)
# print(theta$cholesky)
init_params_j <- matrix(c((theta$D))) # assume a common variance 

# nll_ranef <- function(parms, bi, prior,post){
#   D<-parms
# random_effects_density <- function(u, G, omega, cholesky_parameterization = FALSE, log = FALSE) {
  
#   if (cholesky_parameterization) {
#     # G is provided as Cholesky factor L, so G_actual = L %*% t(L)
#     L <- G
#     G_actual <- (L %*% t(L))
#   } else {
#     # G is the actual covariance matrix
#     G_actual <- G
#   }
  
#   # f(u) = (2π)^(-q/2) |G|^(-1/2) exp(-1/2 u' G^(-1) u)
  
#   if (log) {
#     q <- length(u)
    
#     if (cholesky_parameterization) {
#       # Efficient computation using Cholesky factor
#       # log|G| = 2 * sum(log(diag(L)))
#       log_det_G <- 2 * sum(log(diag(L)))
      
#       # Solve L w = u for w, then u' G^(-1) u = w' w
#       w <- solve(L, u)
#       quadratic_form <- sum(w^2)
      
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
#     return(dnorm(u, mean = rep(0, ncol(u)), sigma = sqrt(G_actual)))
#   }
# }
# cholesky <- chol(D)
# ll <- 0
# # omegas <- c(omegas,1)
# k <- length(bi)
# for(i in 1:nrow(bi[[1]])){
#   for(j in 1:k){
#     ll <- ll+post[i,j]*(log(prior[i,j])+random_effects_density(u=bi[[j]][i],G=t(cholesky), cholesky_parameterization = T, log = T))
#   }
# }
#   return(-ll)
# }

# nll_ranef(D, bi=random_effects, prior=new_prior, post=postprob)

# opt_result <- marqLevAlg::mla(b=D, fn=nll_ranef, bi=bi, prior=new_prior, post=postprob)
# this is increases the log-likelihood (nll should be positive)
# neg_loglik_Dg(init_params_j,data=data,  ranef=random_effects, prior=new_prior, posterior=postprob)


init <- c(theta$betas, log(theta$cholesky),log(theta$sigma_e)) 
nll_marg <- function(parms, data, posterior){
  

marginal_density_yis <- function(data_long, betas_matrix, cholesky, sigmas, D) {
  N <- length(unique(data_long$ID))# number of unique subjects under study.
  m <- nrow(betas_matrix); # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  nis<-table(data_long$ID)
  # if(class(omegas)=="list"){
  #   omegas <- unlist(omegas)
  # }
  # omegas <- c(omegas,1)
  D_matrices <- cholesky
  for (k in 1:m) {
    # k <- 1
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    # if(length(omegas)==1 & k==m){
    #   omega <- 1
    # }else if(length(omegas)==1 & k==1){
    #   omega <- omegas
    # }else if(k<m & length(omegas)!=1){
    # omega <- omegas[k]
    # }else{
    #   omega<-1
    # }# random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects

    for (i in 1:N) {
      # i <- 1
      n_i <- nis[i] # determine number of longitudinal responses per subject
      # k <- 1
      # data_long <- data;i <- 1
      subset <- data_long[data_long$ID == i, ] #! ensure that IDs are consistently named across datasets
      X_i <- model.matrix(~1+measure_time,subset) # Fixed Effects Design Matrix
      p <- ncol(X_i)
      Xi_beta_k <- matrix(X_i %*% beta_g) # Fixed Effects Linear Predictor
      Z_i <- model.matrix(~1, data = subset) # Random Effects Design Matrix
      q <- ncol(Z_i)
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      V_ig <- as.matrix(Z_i %*% (D_g^2) %*% t(Z_i) + (sigma_g^2) * diag(n_i))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig, log = T)
    }
    }
  return(f_yi)
}

  init <- parms
  betas_matrix <- matrix(init[1:4], nrow = 2, ncol=2)
  cholesky <- exp(init[5])
  # cholesky[1,2] <- 0
  sd_e <- exp(init[6])
  # omegas <- init[11:12]
  ysu <- marginal_density_yis(
  data_long = data, 
  betas_matrix = betas_matrix, 
  cholesky= cholesky, 
  sigmas = sd_e)
  
  llk <- sum(rowSums(posterior*(ysu)))
  return(-llk)
}
nll_marg(init, data, postprob)
library(optimParallel)
n_cores <- detectCores() - 1
# opt_result <- marqLevAlg::mla(b=init, fn=nll_marg,data=data, posterior=postprob,nproc=n_cores, print.info=T)
#  opt_result <- obj
cl <- makeCluster(n_cores)     # set the number of processor cores
 setDefaultCluster(cl=cl)
opt_result <- optimParallel(init, fn=nll_marg, data=data, posterior=postprob,hessian=T, parallel = list(cl=cl))
sd_info <- summary(opt_result, loglik=T)
# obj <- opt_result
 (betas_matrix <- matrix(opt_result$par[1:4], nrow = 2, ncol=2,byrow = T))
  cholesky <- exp(metafor::vec2mat(opt_result$par[5], diag=T))
  # cholesky[1,2] <- 0
   (D <- tcrossprod((cholesky)))
  (sd_e <- exp(opt_result$par[6]))
  # (omegas <- opt_result$par[11:12])
cat("\nRE Variance:", D)
cat("\nResidual Variance: ",sd_e^2)
  
# library('optimParallel')
# n_cores <- detectCores() - 1
# # cl <- makeCluster(n_cores)     # set the number of processor cores
# # setDefaultCluster(cl=cl)
# # clusterExport(cl,"dmvnorm", envir = environment())
# # library(tictoc)
# # tic()
# # opt_result <- optimParallel(as.vector(init_params_j), neg_loglik_Dg, ranef=random_effects, prior=new_prior, posterior=postprob,hessian = T, parallel = list(cl=cl, loginfo=T), data=data, verbose=T)
# # toc()
# # stopCluster(cl)

# D <-opt_result$b
# D_g <- matrix(c(parms_j[4], parms_j[5], parms_j[5], parms_j[6]), nrow = T)
# omegas <- opt_result$par[1:(m-1)]
# # print(cholesky)
#     if(!is.null(opt_result$hessian)& !is.singular.matrix(opt_result$hessian)){
# D_sds <- sqrt(diag(matrix.inverse(opt_result$hessian)))

# chol_sd <- metafor::vec2mat(sds[m:(m+2)], diag=T)
# chol_sd[1,2] <- 0
# omegas_sd <- sds[1:(m-1)]
#   } else{
#       omegas_sd <- NULL
#       chol_sd <- NULL
#     }


#   dgs_sds[[k]] <- metafor::vec2mat(sqrt(diag(matrix.inverse(opt_result$hessian))), diag = T)

  # 
# neg_loglik_omegas_Dg(init_params_j,data=data,  ranef=random_effects, prior=new_prior, posterior=postprob)
# cat("Contrib. Omegas:",-neg_loglik_omegas_Dg(matrix(c(sqrt(omegas),cholesky[lower.tri(theta$cholesky, diag = T)])),data=data,  ranef=random_effects, prior=new_prior, posterior=postprob))



#  print("Variance Component Matrices and Associated Standard Errors:")
#  print(D)
# # print(D_sds)
# # print(tcrossprod(cholesky))

# # print("Group specific proportions for Variance-Component Matrices and Associated Standard Errors:")
# # print(sqrt(omegas^2))
# # print(omegas_sd)
# sd_e <- theta$sigma_e; sd_e

# #? compute group specific betas and error variance (this takes an age to compute)
# #' Alternative: Using individual normal densities (more stable)

# #' Total conditional negative log-likelihood (sum over all individuals)
# total_nll_betas <- function(parms, y, random_effects, X, Z, postprob, prior, data) {

# lli_stable <- function(y_i, u_i, X_i, Z_i, beta, sigma_e) {
  
#   # Conditional mean
#   mu_i <- X_i %*% beta + Z_i %*% u_i
#   Sigma <- diag(sigma_e^2, nrow = length(mu_i))
#   # Sum of negative log-densities
#   ll_i <- mvnfast::dmvn(t(y_i), mu = mu_i, sigma = Sigma, log = T)
  
#   return(ll_i)
# }
#   l <- length(parms)
#     betas_matrix <- matrix(parms[-l],ncol = 2) # first p*m parameters are the betas
#   print(betas_matrix)
#   sd_e <- sqrt(parms[l]^2) #the p*m+1 parameter should be the standard error
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
#     X_i <- X[idx_i,, drop=FALSE]
#     Z_i <- Z[idx_i,, drop=FALSE]
#     u_i <- u_matrix[i]  # ith row = [u₀ᵢ, u₁ᵢ]
    
#     # Compute likelihood for individual i
#     ll <- ll + postprob[i,k]*(lli_stable(y_i, u_i = u_i, X_i, Z_i, beta =betas_matrix[k,], sigma_e = sd_e))

#     # nll_i <- lli_stable(y_i, u_i, X_i, Z_i, beta, sigma_e)
#     # total_nll <- total_nll + nll_i
#   }}
#   # total_nll <- sum(rowSums(loglik))
  
#   return(-ll)
# }
# init_params_j <- c(theta$betas, theta$sigma_e)
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

shapes <- as.matrix(theta$shape)
scales <- as.matrix(theta$scale)
gammas_matrix <- theta$gammas
#? compute random effects matrices (D)


# betas_matrix <- matrix(opt_result$par[-5],ncol=2)
# sd_e <- opt_result$par[5]
# # abs(init_pault$par[9]; 
# # ?optim

# # if()
# if(!is.null(opt_result$hessian)& !is.singular.matrix(opt_result$hessian)){
# betas_matrix_sds <- matrix(sqrt(diag(matrix.inverse(opt_result$hessian)))[-5], ncol=2); cat("\nBetas:",betas_matrix,"\nStandard Errors:",betas_matrix_sds)
# sd_e_sd <- sqrt(diag(matrix.inverse(opt_result$hessian)))[5]; cat("\nError Standard Deviation:", sd_e,"\nStandard Errors:",sd_e_sd)
# #! update names
# } else{
   betas_matrix_sds <- NULL
   sd_e_sd <- NULL
# }
# # if(!is.null(opt_result$hessian)& !is.singular.matrix(opt_result$hessian)){
# betas_matrix_sds <- matrix(sqrt(diag(matrix.inverse(opt_result$hessian)))[-7], ncol=2); 
cat("\nBetas:",betas_matrix,"\nStandard Errors:",betas_matrix_sds)
# sd_e_sd <- sqrt(diag(matrix.inverse(opt_result$hessian)))[7];
cat("\nError Standard Deviation:", sd_e,"\nStandard Errors:",sd_e_sd)
# #! update names
# # } else{
#   betas_matrix_sds <- NULL
#   sd_e_sd <- NULL
# # }
# init_params_j <- c(betas_matrix, sd_e)

# cat("\nLog-likelihood Contribution - Betas and Sigma:\n",-total_nll_betas(init_params_j,data=data,postprob=postprob,prior=pi_init, random_effects=random_effects, y=y, X=X, Z=Z))
# # beta_matrix <- beta_gg
# D_matrices <- tcrossprod(cholesky)
# cholesky <- cholesky
sigmas <- sd_e
pi_init <- new_prior

shapes <- as.matrix(theta$shape)
scales <- as.matrix(theta$scale)
# gammas_matrix <- theta$gammas
# wx <- cbind(data_wide$X1, data_wide$X2)
# posterior <- c
# shapes_sd <- shapes
# scales_sd <- scales
# gammas_sd <- gammas_matrix
wb <- density_weibull(shapes, scales)
shapes
scales
# gammas_matrix
parms <- c(log(shapes), log(scales))
neg_loglik_wb <- function(parms_j, data_wide, posterior){
  # parms_j <- parms
  log_shapes <- matrix(parms_j[1:m],ncol = 1); log_shapes # first m parameters are the shapes
  log_scales <- matrix(parms_j[(m+1):(2*m)]); log_scales # next m parameters are the scales
  # gammas_matrix <- matrix(parms_j[(2*m+1):(length(parms_j))], nrow = 4)

  f_td <-  function(shape, scale, t, d) {
    (d*(dweibull(t,shape=shape, scale = scale,log = T))) + (1-d)*(pweibull(t, shape = shape, scale =scale,lower.tail = F,log.p = T)) 
  }

  density_weibull <- function(shape_params, scale_params) {
    N <- nrow(data_wide)
    k <- nrow(as.matrix(shape_params))

    f_td <- matrix(0,nrow = N, ncol = k)
    # w <- model.matrix(~ -1+X3+X4+X5, data=data_wide)
    delta <- data_wide$delta
    # initial density function for Weibull
    for (j in 1:k) {
      # j <- 1
      shape <- shape_params[j]
      scale <- scale_params[j]
      # gammas <- gammas_matrix[j, ]
      delta <- delta
      # fixed effects in survival component
      # wxi <- as.matrix(w)
      ti <- data_wide$event_time
      f_td[, j] <- f_td(shape, scale, ti, delta)
    }
    return(f_td)
  }
  scales <- exp(log_scales)
  shapes <- exp(log_shapes)
  lwb <- density_weibull(shapes, scales)
  component_likelihood <- sum(rowSums(posterior * (lwb)))
  return(-component_likelihood)
}
neg_loglik_wb(parms, data_wide = data_wide, posterior = postprob)
opt_result <- tryCatch({
 optim(par=parms, fn=neg_loglik_wb, data_wide=data_wide, posterior=postprob, hessian=T)}, error = function(e) {
    optim(par=parms, fn=neg_loglik_wb, data_wide=data_wide, posterior=postprob)})


summary(opt_result, loglik=T)
shapes <- exp(matrix(opt_result$par[1:m],ncol = 1)) # first m bameters are the shapes
scales <- exp(matrix(opt_result$par[(m+1):(2*m)])) # next m bameters are the scales
# gammas_matrix <- matrix(opt_result$par[(2*m+1):(length(init_parms_j))], nrow = 4)
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
  # }

cat("\nWeibull Shapes:", shapes, "\nStandard Errors:", shapes_sd)
cat("\nWeibull Scales:", scales, "\nStandard Errors:", scales_sd)

parms <- c(log(shapes), log(scales))
cat("\nContribution for Weibull:", -neg_loglik_wb(parms, data_wide = data_wide, posterior = postprob))
# -neg_loglik_wb(parms, data_wide = data_wide, posterior = postprob)
  
  
bi <- list()
#? compute class specific random effects 
# dgs <- tcrossprod(cholesky)
m <- length(theta$shape)
# omegas_tmp <- c(omegas,1)
for(j in 1:m){
  #  j <- 1
    sigma_g <- sd_e # error variances 

    D_g <- D
     # random effect variance-covariance matrices
    beta_g <- betas_matrix[j,] # fixed effects

  #  i<-1
  bigs<-matrix(nrow = n, ncol = q)
    for(i in 1:n){
    #   i<-2
      n_i <- nis[i]
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


standard_errors <- list()
# update parameters
 theta <- list(pj=new_prior, scale=scales, shape=shapes, betas=betas_matrix, ranef=bi, D=D, sigma_e=sd_e, cholesky=cholesky)
#  standard_errors <- list("Variance Components"=chol_sd, "Omegas"=omegas_sd, "Betas"=betas_matrix_sds, "Error Variance"=sd_e_sd, "Shapes"=shapes_sd, "Scales"=scales_sd)
 ret <- list("Parameters"=theta, "Standard Errors"=standard_errors)
return(ret)
}


# pj <- estep.result$P
# # # # # library(tictoc)
# # # # # tic()
# # # # library(matrixcalc)
#  m_step_result <- mstep(postprob=pj,theta)
# # cat("\n")
# # toc()
# m_step_result$Parameters$pj

# estep.result <- estep(prior_mixing=m_step_result$Parameters$pj, m_step_result$Parameters)
# postprob <- estep.result$P

#  m_step_result <- mstep(postprob=postprob,m_step_result$Parameters)
