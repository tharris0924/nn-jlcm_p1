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
nnem = function(k) {
  # Ensure x is a matrix
  # seed <- 1
  # n <- 250
  # n <<- n
  # seed <<- seed
  # cat("Simulation no.:",seed, "\n")
  # cat("Sample size:",n,"\n")
  # # rm(seed)
  k <- k
  source("NNEM/Application/gen_data.R")
  source("NNEM/Application/em_utils.R")
  source("NNEM/Application/estep.R")
  source("NNEM/Application/mstep.R")
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
    e_step_result <- estep(prior_mixing=pj, theta, long_form_o, long_form_g, ranform, survform, data_long, data_wide)
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
    mstep.result <- mstep(postprob=P,theta, data_wide, data_long, ranform, long_form_g, long_form_o, gform, survform)
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

result <- nnem(k=3)
save(result,file="NNEM/Application/nnjlcmk3.RData")
result <- nnem(k=4)
save(result,file="NNEM/Application/nnjlcmk4.RData")

BIC4=-2*result$loglik+(20*log(2213))

load("NNEM/Application/nnjlcmk2.RData")
ll1 <- result$loglik
(BIC2=(-2*result$loglik)+(20*log(490)))
(AIC2=(-2*result$loglik)+(20*2))

#    AIC: 13331.93  
    #  BIC: 13455.44  

load("NNEM/Application/nnjlcmk3.RData")
k <- nrow(result$theta_sd$`Longitudinal Submodel`)+nrow(result$theta_sd$`Survival Submodel`)
ll2 <- result$loglik
(BIC3=-2*result$loglik+((k)*log(499)))
(AIC3=-2*result$loglik+(25*2))

load("NNEM/Application/nnjlcmk4.RData")
ll4 <- result$loglik
BIC4=-2*result$loglik+(30*log(490))
AIC4=-2*result$loglik+(30*2)

df <- data.frame(ng=c(2,3,4),npm=c(20,25,30),loglik=c(ll1,ll2,ll4),BIC=c(BIC2,BIC3,BIC4), AIC=c(AIC2,AIC3,AIC4))
df
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

