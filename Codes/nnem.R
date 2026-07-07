##Machine l
libs <- c("nnet", "plyr", "metafor", "matrixcalc", "flexsurv", "lcmm", "data.table", "optimParallel", "marqLevAlg", "mvnfast")
install_load <- function (package1, ...)  {   

   # convert arguments to vector
   packages <- c(package1, ...)

   # start loop to determine if each package is installed
   for(package in packages){

       # if package is installed locally, load
       if(package %in% rownames(installed.packages()))
          do.call('library', list(package))

       # if package is not installed locally, download, then load
       else {
          install.packages(package)
          do.call("library", list(package))
       }
   } 
}
install_load(libs)
long_form_g=y~1+ measure_time +(measure_time2)
ranform=~1+ measure_time +(measure_time2)
survform <- Surv(event_times, delta)~1
ranform=~1+ measure_time +(measure_time2)
nnem = function(k, long_form_g, ranform, survform, data_long, data_wide, id_var, y_var, event_var, event_time_var, measure_var, gform) {
  # some small data preprocessing
  data_long$y <- data_long[,y_var]
  data_long$delta <- data_long[,event_var]
  data_long$measure_time <- data_long[,measure_var]
  data_long$measure_time2 <- data_long[,measure_var]^2
  data_long$event_time <- data_long[,event_time_var]
  cat("Initializing parameters...\n")
  # fit a single jlcmm object to get the initial values
  mj1 <- jlcmm(fixed = long_form_g,
               random = ranform,
               survival = Surv(time=event_time, event=delta)~1, hazard = "Weibull",
               subject = id_var, data = data_long, verbose = F, logscale = T, maxiter=250)
  summary(mj1)
  
  
  init <- estimates(mj1, cholesky=F)

  vcovmj1 <- vcov(mj1)
  
  sample_init <- mvnfast::rmvn(k, init, sigma = vcovmj1)
  colnames(sample_init) <- names(init)
  
  D <- sample_init[1,which(grepl('varcov', colnames(sample_init)))]
  cholesky <- metafor::vec2mat(D, diag=T)
  cholesky[3,1] <-cholesky[2,2] 
  cholesky[2,2] <- cholesky[1,3]
  cholesky[1,3] <- cholesky[3,1]
  
  idWeib_shape <- which(grepl('Weibull2', colnames(sample_init)))
  idWeib_scale <- which(grepl('Weibull1', colnames(sample_init)))
  
  (scale_est <- (exp(sample_init[,idWeib_scale]))^(-1/exp(sample_init[,idWeib_shape])))
  (shape_est <- exp(sample_init[,idWeib_shape]))

  
  idxs_slope <- which(grepl(c('intercept', "measure_time", "measure_time2")[2], colnames(sample_init)))
  
  b1 <- sample_init[,(idxs_slope[1]-1)]
  b2 <- sample_init[,idxs_slope[1]]
  b3 <- sample_init[(idxs_slope[2])]
  
  df <- df2 <- data.frame(int=b1,mt=b2,mt2=b3)

  betas_matrix <- as.matrix(df)
  
  residual_err <-sample_init[1,"stderr"]
  
  pi_init <- matrix(rep(1/k), nrow=nrow(data_wide), ncol=k)
  
  cholesky <- t(chol(cholesky))
  # omegas <- alphas[1:2]
 theta <- list(pj=pi_init, scale=scale_est, shape=shape_est, 
              betas=betas_matrix, D=D, cholesky=cholesky, sigma_e=residual_err)
  cat("Initialization complete!\n")
  # Ensure x is a matrix
  # seed <- 1
  # n <- 250
  # n <<- n
  # seed <<- seed
  # cat("Simulation no.:",seed, "\n")
  # cat("Sample size:",n,"\n")
  # # rm(seed)
  k <- k
  # source("NNEM/Application/gen_data.R")
  source("~/Documents/GitHub/nn-jlcm_p1/Application/paquid/em_utils.R")
  source("~/Documents/GitHub/nn-jlcm_p1/Application/paquid/estep.R")
  source("~/Documents/GitHub/nn-jlcm_p1/Application/paquid/mstep.R")
  # print("Initialization complete!")
  # print(theta)
  # theta_orig <- theta
  # theta <- theta_orig
  # Initialize parameters if not provided
  # i can show two instances where the EM algorithm operates but unfortunately
  # does not converge to the "correct" values. Perhaps we should use alternate methods for 
  # estimating the parameters
  # Initialize with provided or computed values
  # m_step <- init
  prev_loglik <- -Inf
  pj <- theta$pj
  # EM iterations
  # iter <- 1
  maxit <- 100
  maxit <- maxit
  # tol <- tol
  # verbose <- T
   tol <- 0.00001
  options(scipen = 999, digits = 6)
  # library(svMisc)
  # pb = txtProgressBar(min = 0, max = maxit, initial = 0)
  # prev_loglik <- e_step_result$loglik
  verbose=T
  options(width = 80)
  likelihoods <- numeric(maxit)
  for (iter in 1:maxit) {
    # print(setTxtProgressBar(pb,iter))
    # progress(iter, progress.bar = T, max.value = maxit)

    # E-step: compute responsibilities (posteriors)
    e_step_result <- estep(prior_mixing=pj, theta = theta, long_form_g = long_form_g, ranform = ranform, data_long = data_long, data_wide = data_wide)
    P <- e_step_result$P
    curr_loglik <- e_step_result$loglik
    likelihoods[iter] <- curr_loglik
    cat("Posterior Probabilities: \n")
    print(round(head(P),digits = 6))
    # if (verbose && iter %% 10 == 0) {
    if (verbose) {
      cat("Iteration", iter, ", log-likelihood:", curr_loglik, "\n")
    }
    
    # Check convergence
    diffs <- abs(curr_loglik - prev_loglik)/abs(prev_loglik); cat("Tolerance check:",diffs,"\n")
    if (iter > 1 && diffs < tol) {
      if (verbose) cat("Converged at iteration", iter, "\n")
      break
    }
    
    prev_loglik <- curr_loglik
    
    
    # M-step: update parameters using the responsibilities
    mstep.result <- mstep(postprob=P,theta = theta, data_wide = data_wide, data_long = data_long, ranform = ranform, long_form_g = long_form_g,gform =  gform)
    theta <- mstep.result$Parameters
    theta_sd <- mstep.result$`Model Estimates & S.E.`
    pj <- theta$pj # update priors

    # progress bar
  #     extra <- nchar('||100%')
  # width <- options()$width
  # step <- round(iter / maxit * (width - extra))
  # text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
  #                 strrep(' ', width - step - extra), round(iter / maxit * 100))
  # cat("\nNNEM Iterations Progress:\n",text, '\n')
  # Sys.sleep(0.05)
  # cat(if (iter == maxit) '\n' else '\014')

  }
  
  # Get final cluster assignments
  C <- apply(P, 1, which.max)
  table(C)
  # # iter_name <- 
  #  filename <- paste0("nnem_size", seed, n,".RData")
  # save(result, file = filename)
  # save(result,)
#   if(verbose){
    
#     # theta <- result$theta
# cat("\nError Variance:\n","Est",theta$sigma_e,theta_orig$sigma_e,(abs(theta$sigma_e - theta_orig$sigma_e)))
# cat("Cholesky");theta_orig$cholesky; theta$cholesky; (bias_cholesky <- abs(theta$cholesky - theta_orig$cholesky))
# cat("RE Variance-Covariance");vech(theta_orig$D);vech(theta$D);(bias_D <- abs(vech(theta$D)-vech(theta_orig$D)))
# cat("Omegas:");sqrt((theta$omegas)^2); theta_orig$omegas;(bias_omegas <- abs(sqrt((theta$omegas)^2) - theta_orig$omegas))
# cat("Longitudinal Fixef");theta$betas; theta_orig$betas;(bias_beta <- abs(theta$betas - theta_orig$betas))
# cat("Weibull Scales"); theta$scale; theta_orig$scale;(bias_scale <- abs((theta$scale)-(theta_orig$scale)))
# cat("Weibull Shapes"); theta$shape; theta_orig$shape; (bias_shape<- abs(theta$shape - theta_orig$shape))

#   }
  # Return results
  k <- nrow(theta_sd$`Longitudinal Submodel`) + nrow(theta_sd$`Survival Submodel`)
  n <- nrow(data_wide)
  BIC <- -2*(curr_loglik)+ k*log(n)
  AIC=(-2*curr_loglik)+(k*2)
  result <- list(
    theta = theta,
    theta_sd = theta_sd,
    pj = theta$pj,
    iteration = iter,
    BIC=BIC,
    AIC=AIC,
    C = C,
    loglik = curr_loglik,
    all_loglik = likelihoods,
    # x=x,
    postprob=P,
    data_wide=data_wide,
    data_long=data
  )
  return(result)
}


# run an instance of the simulation from the paper K=4
source("gen_data_k4.R")
result <- nnem(k=4, long_form_g = long_form_g, ranform = ranform, data_long = data_long, data_wide = data_wide, id_var = "ID", gform=~x, y_var = "y", event_var = "delta", measure_var = "measure_time", event_time_var = "event_time")

# run the model on the application dataset, PAQUID
source("gen_data_paquid.R")
result <- nnem(k=4, long_form_g = long_form_g, ranform = ranform, data_long = data_long, data_wide = data_wide, id_var = "ID", gform=~x, y_var = "y", event_var = "delta", measure_var = "measure_time", event_time_var = "event_time")
