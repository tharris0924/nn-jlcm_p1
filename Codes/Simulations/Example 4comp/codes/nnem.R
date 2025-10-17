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
nnem = function(n, seed, init = NULL, maxit = 100, tol = 1e-4, verbose = TRUE) {
  # Ensure x is a matrix
  # seed <- 1
  # n <- 250
  n <<- n
  seed <<- seed
  cat("Simulation no.:",seed, "\n")
  cat("Sample size:",n,"\n")
  # rm(seed)
  source("NNEM/Example_k4/gen_data.R")
  source("NNEM/Example_k4/em_utils.R")
  source("NNEM/Example_k4/estep.R")
  source("NNEM/Example_k4/mstep.R")
  print("Initialization complete!")
  # print(theta)
  theta_orig <- theta
  # Initialize parameters if not provided
  # i can show two instances where the EM algorithm operates but unfortunately
  # does not converge to the "correct" values. Perhaps we should use alternate methods for 
  # estimating the parameters
  # Initialize with provided or computed values
  # m_step <- init
  prev_loglik <- -Inf
  pj <- pi_init
  # EM iterations
  # iter <- 1
  maxit <- 100
  maxit <- maxit
  tol <- tol
  # verbose <- T
  # tol <- 0.001
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
    e_step_result <- estep(pj, theta)
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
    mstep.result <- mstep(P, theta)
    theta <- mstep.result$Parameters
    theta_sd <- mstep.result$`Standard Errors`
    pj <- theta$pj # update priors

    # progress bar
      extra <- nchar('||100%')
  width <- options()$width
  step <- round(iter / maxit * (width - extra))
  text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
                  strrep(' ', width - step - extra), round(iter / maxit * 100))
  cat("\nNNEM Iterations Progress:\n",text, '\n')
  Sys.sleep(0.05)
  cat(if (iter == maxit) '\n' else '\014')

  }
  
  # Get final cluster assignments
  C <- apply(P, 1, which.max)
  table(C)
  # # iter_name <- 
  #  filename <- paste0("nnem_size", seed, n,".RData")
  # save(result, file = filename)
  # save(result,)
  if(verbose){
    
    # theta <- result$theta
cat("\nError Variance:\n","Est",theta$sigma_e,theta_orig$sigma_e,(abs(theta$sigma_e - theta_orig$sigma_e)))
cat("Cholesky");theta_orig$cholesky; theta$cholesky; (bias_cholesky <- abs(theta$cholesky - theta_orig$cholesky))
cat("RE Variance-Covariance");vech(theta_orig$D);vech(theta$D);(bias_D <- abs(vech(theta$D)-vech(theta_orig$D)))
cat("Omegas:");sqrt((theta$omegas)^2); theta_orig$omegas;(bias_omegas <- abs(sqrt((theta$omegas)^2) - theta_orig$omegas))
cat("Longitudinal Fixef");theta$betas; theta_orig$betas;(bias_beta <- abs(theta$betas - theta_orig$betas))
cat("Weibull Scales"); theta$scale; theta_orig$scale;(bias_scale <- abs((theta$scale)-(theta_orig$scale)))
cat("Weibull Shapes"); theta$shape; theta_orig$shape; (bias_shape<- abs(theta$shape - theta_orig$shape))

  }
  # Return results
  result <- list(
    theta = theta,
    theta_sd = theta_sd,
    pj = theta$pj,
    iteration = iter,
    C = C,
    loglik = curr_loglik,
    all_loglik = likelihoods,
    x=x,
    data_wide=data_wide,
    data_long=data
  )
  return(result)
}
# n <- 1000
# 0.880044 1 0.119956

# n <- 500
# 0.882017 1 0.117983

# n <- 100
#  0.867987 1 0.132013

#  Est 0.954322 1 0.0456783
# 0.975474 1 0.0245263

# Est 0.873613 1 0.126387
#  Est 0.873613 1 0.126387
n <- 100
seed <- 8
censor_rate <- 0.05 
  # r <- 0.05 # 5%
  # # r <- 0.2 # 25%
  # # r <- 0.5 # 50%
# rm(list=setdiff(ls(), c("seed","n", "censor_rate", "nnem","wkdir")))
sdsds <- nnem(n=n,seed, tol = 0.0001 ) #0.92071

M <- 500
# we need to run for up to 5000 subjects and we need to demonstrate model performance under 0% censoring
for(cens in c(0.05, 0.25, 0.5)){
censor_rate <- cens
names <- "Results_n=250"
wkdir <- paste0("k=3_",names,"_c",censor_rate)
dir.create(wkdir)
for(seed in 1:(M)){
  # burn everything
  print(seed)
  # rm(list=ls())
  rm(list=setdiff(ls(), c("seed","n", "nnem","wkdir", "cens", "censor_rate", "M")))
  # seed <- 116
  result <- tryCatch({nnem(n, seed, tol = 0.0001)}, error=function(e){NULL})
  source("NNEM/Example_k4/jlcmm_fit.R")
  filename <- paste0(wkdir,"/nnem_size", seed, n,".RData")
  save(result, file = filename)
  filename2 <- paste0(wkdir,"/jlcmm_fit", seed, n,".RData")
  save(jlcmm_fit, file = filename2)
}
}

n <- 1000
# seed <- 8
# censor_rate <- 0.05 
  # r <- 0.05 # 5%
  # # r <- 0.2 # 25%
  # # r <- 0.5 # 50%
# rm(list=setdiff(ls(), c("seed","n", "censor_rate", "nnem","wkdir")))
M <- 500
# we need to run for up to 5000 subjects and we need to demonstrate model performance under 0% censoring
for(cens in c(0.05, 0.25, 0.5)){
censor_rate <- cens
# sdsds <- nnem(n=n,seed, tol = 0.000001 ) #0.92071
names <- "Results_n=1000"
wkdir <- paste0("k=3_",names,"_c",censor_rate)
dir.create(wkdir)
for(seed in 1:(M)){
  # burn everything
  print(seed)
  # rm(list=ls())
  rm(list=setdiff(ls(), c("seed","n", "nnem","wkdir", "cens", "censor_rate", "M")))
  # seed <- 116
  result <- tryCatch({nnem(n, seed, tol = 0.0001)}, error=function(e){NULL})
  source("NNEM/Example_k-3/jlcmm_fit.R")
  filename <- paste0(wkdir,"/nnem_size", seed, n,".RData")
  save(result, file = filename)
  filename2 <- paste0(wkdir,"/jlcmm_fit", seed, n,".RData")
  save(jlcmm_fit, file = filename2)
}
}

# load("nnem_size250.RData")
theta <- result$theta
print("Error Variance:");theta$sigma_e;theta_orig$sigma_e;(bias_sigma <- abs(theta$sigma_e - theta_orig$sigma_e))
print("Cholesky");theta_orig$cholesky; theta$cholesky; (bias_cholesky <- abs(theta$cholesky - theta_orig$cholesky))
print("RE Variance-Covariance");vech(theta_orig$D);vech(theta$D);(bias_D <- abs(vech(theta$D)-vech(theta_orig$D)))
print("Omegas:");sqrt((theta$omegas)^2); theta_orig$omegas;(bias_omegas <- abs(sqrt((theta$omegas)^2) - theta_orig$omegas))
print("Longitudinal Fixef");theta$betas; theta_orig$betas;(bias_beta <- abs(theta$betas - theta_orig$betas))
print("Weibull Scales"); theta$scale; theta_orig$scale;(bias_scale <- abs((theta$scale)-(theta_orig$scale)))
print("Weibull Shapes"); theta$shape; theta_orig$shape; (bias_shape<- abs(theta$shape - theta_orig$shape))

