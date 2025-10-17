# density of survival times with censoring
  f_td <-  function(shape, scale, t, d) {
    ((dweibull(t,shape=shape, scale = scale)^d)) * (pweibull(t, shape = shape, scale =scale,lower.tail = F)^((1-d))) 
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


marginal_density_yis <- function(data_long, betas_matrix, cholesky, sigmas, omegas) {
  N <- length(unique(data_long$ID));print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix);print(m);print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  nis<-table(data_long$ID)
  if(class(omegas)=="list"){
    omegas <- unlist(omegas)
  }
  omegas <- c(omegas,1)
  D_matrices <- tcrossprod(cholesky)
  for (k in 1:m) {
    # k <- 1
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    # if(length(omegas)==1 & k==m){
    #   omega <- 1
    # }else if(length(omegas)==1 & k==1){
    #   omega <- omegas
    # }else if(k<m & length(omegas)!=1){
    omega <- omegas[k]
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
      Z_i <- X_i # Random Effects Design Matrix
      q <- ncol(Z_i)
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      V_ig <- as.matrix(Z_i %*% ((omega^2)*D_g) %*% t(Z_i) + (sigma_g^2 * diag(n_i)))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig)
    }
    }
  return(f_yi)
}




# # example usage
# ysu <- marginal_density_yis(
#   data_long = data, 
#   betas_matrix = betas_matrix, 
#   cholesky= cholesky, 
#   sigmas = sd_e, 
#   omegas = theta$omegas)

# log_sde <- log(1)
# init <- c(betas_matrix, cholesky[lower.tri(cholesky, diag = T)],log_sde, omegas) 
# nll_marg <- function(parms, data, posterior){
  

# marginal_density_yis <- function(data_long, betas_matrix, cholesky, sigmas, omegas) {
#   N <- length(unique(data_long$ID)); # number of unique subjects under study.
#   m <- nrow(betas_matrix); # number of unique groups for mixture model
#   f_yi <- matrix(0,nrow = N, ncol = m)
#   nis<-table(data_long$ID)
#   if(class(omegas)=="list"){
#     omegas <- unlist(omegas)
#   }
#   omegas <- c(omegas,1)
#   D_matrices <- tcrossprod(cholesky)
#   for (k in 1:m) {
#     # k <- 1
#     sigma_g <- sigmas # error variances 
#     D_g <- D_matrices # random effect variance-covariance matrices
#     # if(length(omegas)==1 & k==m){
#     #   omega <- 1
#     # }else if(length(omegas)==1 & k==1){
#     #   omega <- omegas
#     # }else if(k<m & length(omegas)!=1){
#     omega <- omegas[k]
#     # }else{
#     #   omega<-1
#     # }# random effect variance-covariance matrices
#     beta_g <- betas_matrix[k,] # fixed effects

#     for (i in 1:N) {
#       # i <- 1
#       n_i <- nis[i] # determine number of longitudinal responses per subject
#       # k <- 1
#       # data_long <- data;i <- 1
#       subset <- data_long[data_long$ID == i, ] #! ensure that IDs are consistently named across datasets
#       X_i <- model.matrix(~1+measure_time,subset) # Fixed Effects Design Matrix
#       p <- ncol(X_i)
#       Xi_beta_k <- matrix(X_i %*% beta_g) # Fixed Effects Linear Predictor
#       Z_i <- X_i # Random Effects Design Matrix
#       q <- ncol(Z_i)
#       Y_i <- subset$y # Longitudinal recordings

#       #? calculating the distributional parameters for the group specific random effects
#       # D_g_inv <- solve(D_g)
#       # 3-by-2 dimensions
#       V_ig <- as.matrix(Z_i %*% ((omega^2)*D_g) %*% t(Z_i) + (sigma_g^2 * diag(n_i)))
#       # simgas[i] <- list(Sigma_i)
#       f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig)
#     }
#     }
#   return(f_yi)
# }
#   init <- parms
#   betas_matrix <- matrix(init[1:6], nrow = 3, ncol=2,byrow = T)
#   cholesky <- metafor::vec2mat(init[7:9], diag=T)
#   cholesky[1,2] <- 0
#   sd_e <- exp(init[10])
#   omegas <- init[11:12]
#   ysu <- marginal_density_yis(
#   data_long = data, 
#   betas_matrix = betas_matrix, 
#   cholesky= cholesky, 
#   sigmas = sd_e, 
#   omegas = omegas)
  
#   llk <- sum(rowSums(posterior*log(ysu)))
#   return(-llk)
# }
# nll_marg(init, data, postprob)

# obj <- mla(init, fn=nll_marg,data=data, posterior=postprob,nproc=19, print.info = T)
# summary(obj, loglik=T)

#  ( betas_matrix <- matrix(obj$b[1:6], nrow = 3, ncol=2,byrow = T))
#   cholesky <- metafor::vec2mat(obj$b[7:9], diag=T)
#   cholesky[1,2] <- 0
#    (D <- tcrossprod(cholesky))
#   (sd_e <- exp(obj$b[10]))
#   (omegas <- obj$b[11:12])
log_conditional_density_yis <- function(data_long, betas_matrix, sigmas, random_effects) {
  #! use this density function to find the betas and error variance
  N <- length(unique(data_long$ID)) #;print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix)#;print(m);print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  for (k in 1:m) {
    # k <- 1
    # random_effects <- bi
    # D_matrices <- D_g_og
    b_g <- as.matrix(random_effects[[k]]) # random effects matrices
    # sigma_g <- sigmas[[k]] # error variances 
    sigma_g <- sigmas # error variances 

    # D_g <- D_matrices[[k]]
     # random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects

    for (i in 1:N) {
      n_i <- nis[i] # determine number of longitudinal responses per subject
      # k <- 1
      # data_long <- sim
      subset <- data_long[data_long$ID == i, ] #! ensure that IDs are consistently named across datasets
      X_i <- model.matrix(~1+measure_time,subset) # Fixed Effects Design Matrix
      p <- ncol(X_i)
      Xi_beta_k <- X_i %*% beta_g # Fixed Effects Linear Predictor
      Z_i <- X_i # Random Effects Design Matrix
      q <- ncol(Z_i)
      eta_i <- Xi_beta_k+Z_i%*%b_g[i,]
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      # Sigma_i <- as.matrix(Z_i %*% D_g %*% t(Z_i) + sigma_g * diag(n_i))
      R_i <- diag(n_i)*sigma_g^2
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- dmvn(t(Y_i), eta_i, R_i, log = T)
    }
  }
  return(f_yi)
}


log_density_weibull <- function(shape_params, scale_params) {
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
    l_i <- (delta*(log(scale*shape)+shape*log(ti))-(scale*ti^shape))
    f_td[, j] <- l_i
  }
  return(f_td)
}

# data_long <- data
# random_effects <- theta$ranef
# sigmas <- sd_e
conditional_density_yis <- function(data_long, betas_matrix, sigmas, random_effects) {
  #! use this density function to find the betas and error variance
  N <- length(unique(data_long$ID)) #;print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix)#;print(m);print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
    nis<-table(data_long$ID)
  # if(class(omegas)=="list"){
  #   omegas <- unlist(omegas)
  # }
  # omegas <- c(omegas,1)
  # D_matrices <- tcrossprod(cholesky)
  for (k in 1:m) {
    # k <- 1
    # random_effects <- bi
    # D_matrices <- D_g_og
    b_g <- as.matrix(random_effects[[k]]) # random effects matrices
    # sigma_g <- sigmas[[k]] # error variances 
    sigma_g <- sigmas # error variances 

    # D_g <- D_matrices[[k]]
     # random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects

    for (i in 1:N) {
      n_i <- nis[i] # determine number of longitudinal responses per subject
      #  i <- 1
      # data_long <- sim
      subset <- data_long[data_long$ID == i, ] #! ensure that IDs are consistently named across datasets
      X_i <- model.matrix(~1+measure_time,subset) # Fixed Effects Design Matrix
      p <- ncol(X_i)
      Xi_beta_k <- X_i %*% beta_g # Fixed Effects Linear Predictor
      Z_i <- X_i # Random Effects Design Matrix
      q <- ncol(Z_i)
      eta_i <- Xi_beta_k+Z_i%*%b_g[i,]
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      # Sigma_i <- as.matrix(Z_i %*% D_g %*% t(Z_i) + sigma_g * diag(n_i))
      R_i <- diag(n_i)*sigma_g^2
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), eta_i, R_i)
    }
  }
  return(f_yi)
}

# data_long <- data

# cond <- conditional_density_yis(data, betas_matrix, sigmas = sd_e, theta$ranef)

density_bis <- function(data_long, cholesky, random_effects, omegas) {
  #! use this density function to find the betas and error variance
  N <- length(unique(data_long$ID)); # number of unique subjects under study.
  m <- length(random_effects) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  if(class(omegas)=="list"){
    omegas <- unlist(omegas)
  }
  # print(cholesky)
  omegas <- c(omegas,1)

  D_matrices <- crossprod(cholesky)
  for (k in 1:m) {
    # k <- 1
    # random_effects <- bi
    # D_matrices <- D_g_og
    b_g <- random_effects[[k]] # random effects matrices
    # sigma_g <- sigmas[[k]] # error variances 
    # sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    # random effect variance-covariance matrices
     omega <- omegas[k]
    # beta_g <- betas_matrix[k,] # fixed effects

    for (i in 1:N) {
      mu_g <- c(0,0)
      f_yi[i, k] <- (mvnfast::dmvn(b_g[i,], mu_g, (omega^2)*D_g))
    }
    }
  return(f_yi)
}

# fd <- density_bis(data, cholesky, theta$ranef, theta$omegas)*conditional_density_yis(data, betas_matrix, sigmas = sd_e, theta$ranef)
