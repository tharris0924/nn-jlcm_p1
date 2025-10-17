# =========================================
# E-step
#=========================================
longform=y~1+ measure_time +(measure_time2) + CEP
ranform=~1+ measure_time +(measure_time2)
data_long <- datalong
survform=~CEP+male
# prior_mixing=result$theta$pj
estep <- function(prior_mixing, theta, long_form_o, long_form_g, ranform, survform, data_long, data_wide) {


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
long_form_o = long_form_o,
ranform = ranform)

  # print(f_yis)
  # f_yis <- density_bis(data, cholesky, theta$ranef, theta$omegas)*conditional_density_yis(data, betas_matrix, sigmas = sd_e, theta$ranef)


  gammas_params <- theta$gamma_matrix
  scale_params <- as.matrix(theta$scale)
  shape_params <- as.matrix(theta$shape)
  # density for time-to-event component
  # f_td <- density_weibull(shapes, scales, gammas_matrix)
    f_td <-  function(shape, scale, t, d) {
    ((dweibull(t,shape=shape, scale = scale)^d)) * (pweibull(t, shape = shape, scale =scale,lower.tail = F)^((1-d))) 
  }

  f_tds <- density_weibull(shape_params = shape_params, scale_params = scale_params, gammas_params = gammas_params, data_wide=data_wide, survform)
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
  # LogLik1=sum(log(rowSums(sapply(1:m,function(j) pj[,j]*f_yis[,j]*f_tds[,j])))); print(LogLik1)

  # Normalize responsibilities
  row_sums <- rowSums(comp)
  comp.p <- comp / row_sums
  
  # Handle any NaN values
  comp.p[is.nan(comp.p)] <- 1/m
  
  # Calculate log-likelihood
  loglik <- sum(log(row_sums)); #print(loglik)
  
  return(list("loglik" = loglik, "P" = comp.p))
}

 
# E-step working!!
# estep.result <- estep(prior_mixing=theta$pj,theta, long_form_o, long_form_g, ranform, survform, data_long, data_wide)

# postprob <- estep.result$P
# loglik <- estep.result$loglik
