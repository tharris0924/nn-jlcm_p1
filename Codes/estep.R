# =========================================
# E-step
#=========================================

# prior_mixing=result$theta$pj
estep <- function(prior_mixing, theta, long_form_g, ranform, data_long, data_wide) {


  shapes <- as.matrix(theta$shape)
  betas_matrix <- theta$betas
  cholesky <- theta$cholesky
  sigmas <- theta$sigma_e

  
  #number of components
  m <- length(shapes)
# density for longitudinal component
  f_yis <- marginal_density_yis(
  data_long = data_long, 
theta=theta,
long_form_g = long_form_g,
ranform = ranform)


  scale_params <- as.matrix(theta$scale)
  shape_params <- as.matrix(theta$shape)
  # density for time-to-event component
    f_td <-  function(shape, scale, t, d) {
    ((dweibull(t,shape=shape, scale = scale)^d)) * (pweibull(t, shape = shape, scale =scale,lower.tail = F)^((1-d))) 
  }

  f_tds <- density_weibull(shape_params = shape_params, scale_params = scale_params, data_wide=data_wide)
  # print(f_tds)
  #
  comp <- matrix(ncol = m, nrow = nrow(data_wide))
  # prior_mixing <- pi_init
  pj <- prior_mixing
  # pj <- prob
  # Calculate component likelihoods
  for (j in 1:m) {
    comp[, j] <-  pj[,j]*f_yis[,j]*f_tds[,j]
  }

  # Normalize responsibilities
  row_sums <- rowSums(comp)
  comp.p <- comp / row_sums
  
  # Handle any NaN values
  comp.p[is.nan(comp.p)] <- 1/m
  
  # Calculate log-likelihood
  loglik <- sum(log(row_sums)); #print(loglik)
  
  return(list("loglik" = loglik, "P" = comp.p))
}

