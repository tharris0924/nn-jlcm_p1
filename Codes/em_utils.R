# density of survival times with censoring
 f_td <-  function(shape, scale, t, d) {
    ((dweibull(t,shape=shape, scale = scale)^d)) * (pweibull(t, shape = shape, scale =scale,lower.tail = F)^((1-d))) 
  }

density_weibull <- function(shape_params, scale_params, data_wide) {

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
    # linear_pred <- exp()
    delta <- data_wide$delta
    # fixed effects in survival component
    # wxi <- as.matrix(w)
    ti <- data_wide$event_time
    for(i in 1:N){
    f_tds[i, j] <- f_td(shape, scale, ti[i], delta[i])
    }
    }
  return(f_tds)
}

# ftds <- density_weibull(shape_params = theta$shape, scale_params = theta$scale, gammas_params = gammas_params, data_wide, ~CEP+male)

# data_long <- datalong
marginal_density_yis <- function(data_long, theta, long_form_g,ranform) {
  # parms prep
  betas_matrix <- theta$betas
  sigmas <- theta$sigma_e
  # D_matrices <- theta$D

  # data prep
  y<- model.frame(long_form_g,data_long)[,1]
  Xg <- model.matrix(long_form_g, data_long)
  Zg <- model.matrix(ranform, data_long)
  data_long$y <- y

  N <- length(unique(data_long$ID)) #;print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix) #;print(m);print(betas_matrix) # number of unique groups for mixture model
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
      Z_ig <- Zg[data_long$ID == ids, , drop = FALSE]    # Random Effects Design Matrix
      X_ig <- Xg[data_long$ID==ids,,drop=FALSE ]# Cluster-specific Fixed Effects Design Matrix

      p <- ncol(X_ig)
      Xi_beta_k <- matrix(X_ig %*% beta_g) # Fixed Effects Linear Predictor
   
      q <- ncol(Z_ig)
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      V_ig <- as.matrix(Z_ig %*% (D_g) %*% t(Z_ig) + sigma_g^2 * diag(nrow(X_ig)))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig)
    }
    }
  return(f_yi)
}
