# =========================================
# E-step
#=========================================
estep <- function(prior_mixing, theta) {


  shapes <- as.matrix(theta$shape)
  betas_matrix <- theta$betas
  cholesky <- theta$cholesky
  sigmas <- theta$sigma_e

  
  #number of components
  m <- length(shapes)
# density for longitudinal component
  f_yis <- marginal_density_yis(
  data_long = data, 
  betas_matrix = betas_matrix, 
  cholesky = cholesky,
  sigmas = theta$sigma_e, 
  omegas = theta$omegas)

  gammas_matrix <- theta$gammas
  scale_params <- as.matrix(theta$scale)
  # density for time-to-event component
  # f_td <- density_weibull(shapes, scales, gammas_matrix)
  f_td <- density_weibull(shapes, scale_params)

  #
  comp <- matrix(ncol = m, nrow = nrow(data_wide))
  # prior_mixing <- pi_init
  pj <- prior_mixing
  # Calculate component likelihoods
  for (j in 1:m) {
    comp[, j] <-  pj[,j]*f_yis[,j]*f_td[,j]
  }
  LogLik1=sum(log(rowSums(sapply(1:m,function(j) pj[,j]*f_yis[,j]*f_td[,j])))); print(LogLik1)

  # Normalize responsibilities
  row_sums <- rowSums(comp)
  comp.p <- comp / row_sums
  
  # Handle any NaN values
  comp.p[is.nan(comp.p)] <- 1/m
  
  # Calculate log-likelihood
  loglik <- sum(log(row_sums)); print(loglik)
  
  return(list("loglik" = loglik, "P" = comp.p))
}


 
# E-step working!!
# estep.result <- estep(prior_mixing=pi_init, theta)
# postprob <- estep.result$P
# # loglik <- estep.result$loglik
