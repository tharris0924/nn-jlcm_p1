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


marginal_density_yis <- function(data_long, betas_matrix, cholesky, sigmas) {
  N <- length(unique(data_long$ID));print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix);print(m);print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  nis<-table(data_long$ID)

  D_matrices <- tcrossprod(cholesky)
  for (k in 1:m) {
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects

    for (i in 1:N) {
      n_i <- nis[i] # determine number of longitudinal responses per subject
      # k <- 1
      # data_long <- data;i <- 1
      subset <- data_long[data_long$ID == i, ] #! ensure that IDs are consistently named across datasets
            # subset<-sim[sim$ID==i,]
      Z_i <- Z[data_long$ID == i, , drop = FALSE]
      X_i <- Z_i# Fixed Effects Design Matrix
      p <- ncol(X_i)
      Xi_beta_k <- matrix(X_i %*% beta_g) # Fixed Effects Linear Predictor
      Z_i <- X_i # Random Effects Design Matrix
      q <- ncol(Z_i)
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      V_ig <- as.matrix(Z_i %*% (D_g) %*% t(Z_i) + sigma_g^2 * diag(n_i))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- dmvn(t(Y_i), Xi_beta_k, V_ig)
    }
    }
  return(f_yi)
}


# example usage
# ysu <- marginal_density_yis(
#   data_long = data, 
#   betas_matrix = betas_matrix, 
#   cholesky= cholesky, 
#   sigmas = sd_e)


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

