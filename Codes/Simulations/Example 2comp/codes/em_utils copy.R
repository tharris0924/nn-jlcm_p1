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
#  ISE_out[i] <- get_nnem_ISE(result, shapes= shape_est[i,], scales=scale_est[i,], mix.prob = comp.p,theta_orig=result$theta_orig, true_pi=true_comp)
# mix.prob <- true_pi
get_nnem_ISE <- function(result, shapes, scales, mix.prob,theta_orig, true_pi){
  ISE <- 0
      times <- seq(0,max(result$data_wide$event_time), length.out = 1000)
    ntimes <- length(times)
      for (x in c(1:nrow(result$data_wide))){
        Shat <- pred_surv_marg(shape_params = shapes, scale_params=scales, mix.prob=t(as.matrix(mix.prob[x,])), data_wide=result$data_wide,x)
        Strue <- pred_surv_marg(shape_params = theta_orig$shape, scale_params=theta_orig$scale,mix.prob=t(as.matrix(true_pi[x,])), data_wide=result$data_wide)
        scores <- (Shat[,2] - Strue[,2])^2
        SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(times)) / diff(range(times))
        ISE <- ISE + SE
    }
    ISE <- ISE / nrow(result$data_wide)
  return(ISE)
}
get_nnem_brier_score <- function(shape_params, scale_params, test_data, mix.prob) {
 #* determine class belongings
  test_data <- data.table(test_data)
  test_data$rownames <- c(1:nrow(test_data))
 rowtokeep <- test_data[, rownames[1], by = ID][, V1]
   test_data <- test_data[rowtokeep,]
  colnames(mix.prob) <- c("V1","V2")
    # pred <- matrix(0, nrow = N, ncol = k)
  library(tidyverse)
  # mix.prob <- pij[[1]]
  dat <- cbind(test_data,mix.prob)
   dat <- dat %>% arrange((event_time))
   mix.prob <- dat |> select(c(V1,V2)) |> as.matrix()
  #  ti <- seq(min(test_data$event_time), max(test_data$event_time), length.out = 1000)
   ti <- c(dat$event_time)
  # mbvars
  # classmb <- paste0("~",mbvars, sep="")
  
      N <- nrow(test_data)
  k <- length(as.matrix(shape_params))

  # initial density function for Weibull
  pred <- matrix(0,nrow=length(ti), ncol=k )

      # ti <- test_data$event_time
      # status <- result$data_wide$delta
  for (j in 1:k) {
    # j <- 1
    shape <- shape_params[j]
    scale <- scale_params[j]

    pred[,j] <- pweibull(ti,shape, scale, lower.tail = F)
  }

 nclasses <- length(shape_params)
#  ti <- seq(0, max(result$data_wide$event_time), length.out = 100)
#  status <- result$data_wide$delta
est_surv <- pred
 #* determine survival curves
 times <- ti
 ntimes <- length(times)
 # ISE <- 0
 predclass <- mix.prob
 if (length(dim(predclass)) == 0) {
   avg <- FALSE
 } else {
   avg <- TRUE
 }
 data <- test_data
 Shat <- NULL
 for (x in c(1:nrow(test_data))) {
    # x <- 1
   if (!avg) {
     tmpc <- predclass[x]
     tmpebx <- exp(sum(pred_slopes[tmpc, ] * data[x, c("X3", "X4", "X5")]))
     Shat <- exp(-tmpebx * mod$predSurv[, paste0("event1.CumRiskFct", tmpc)])
   } else {
     Shat_raw <- matrix(0, ncol = nclasses, nrow = ntimes)
          # tmpc <- predclass[x]

     for (tmpc in c(1:nclasses)) {
        # tmpc <- 1
      #  tmpebx <- exp(sum(pred_slopes[tmpc, ] * data[x, ..survvars])) 
       Shat_raw[, tmpc] <- est_surv[, tmpc] #* multiplying with the baseline cumulative hazard function estimated from the model:
       #* NOTE, this is the same for all i. What makes this change is when it is multiplied with the mixing proportions below
     }
     Shat <- rbind(Shat, c(Shat_raw %*% predclass[x, ]))
   }
 }

 #* calculate integrated Brier score
#  evaltimes <- seq(min(data$event_time), max(data$event_time), length.out = 1000)
#  torm <- c(which(evaltimes < min(test_data$agedem)), which(evaltimes > max(test_data$agedem)))
 test_data$rownames <- c(1:nrow(test_data))
#  rowtokeep <- test_data[, rownames[1], by = ID][, V1]
#  SHAT_k <- Shat[rowtokeep,]
#  test_data <- test_data[rowtokeep,]
#  SHAT_k <- SHAT_k[, -torm]
 brier_surv <- Surv(test_data$event_time, test_data$delta)
 IBS <- sbrier(obj = brier_surv,pred =  t(Shat), btime =  ti)
  return(list(IBS=IBS,Shat=Shat))
}


marginal_density_yis <- function(data_long, betas_matrix, cholesky, sigmas, D) {
  N <- length(unique(data_long$ID));print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix);print(m);print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  nis<-table(data_long$ID)
  # if(class(omegas)=="list"){
  #   omegas <- unlist(omegas)
  # }
  # omegas <- c(omegas,1)
  D_matrices <- D
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
      V_ig <- as.matrix(Z_i %*% (D_g) %*% t(Z_i) + (sigma_g^2) * diag(n_i))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig)
    }
    }
  return(f_yi)
}


# example usage
# ysu <- marginal_density_yis(
#   data_long = data, 
#   betas_matrix = betas_matrix, 
#   cholesky= cholesky, 
#   sigmas = sd_e, 
#   D=D_g)


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
      Z_i <- model.matrix(~1, data = subset) # Random Effects Design Matrix
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


#* I believe that this is a reimplementation from the ipred package
sbrier <- function(obj, pred, btime = range(obj[,1])){
   if(!inherits(obj, "Surv"))
       stop("obj is not of class Surv")

   # check for right censoring

   # <FIXME> why
   class(obj) <- NULL
   # </FIXME>
   if (attr(obj, "type") != "right")
       stop("only right-censoring allowed")
   N <- nrow(obj)

   # get the times and censoring of the data, order them with resp. to time

  time <- obj[,1]
  ot <- order(time)
  cens <- obj[,2]
  time <- time[ot]

   # get the times to compute the (integrated) Brier score over

   if (is.null(btime)) stop("btime not given")
   if (length(btime) < 1) stop("btime not given")

   if (length(btime) == 2) {
       if (btime[1] < min(time)) warning("btime[1] is smaller than min(time)")
       if (btime[2] > max(time)) warning("btime[2] is larger than max(time)")
       btime <- time[time >= btime[1] & time <=
                                         btime[2]]
   }

   ptype <- class(pred)
   # <begin> S3 workaround
   if (is.null(ptype)) {
     if (is.vector(pred)) ptype <- "vector"
     if (is.list(pred)) ptype <- "list"
   }
   # <end>
   # this code breaks the whole thing so lets remove it
  #  if (ptype == "numeric" && is.vector(pred)) ptype <- "vector"

   survs <- NULL
   switch(ptype[1], survfit = {
       survs <- getsurv(pred, btime)
       survs <- matrix(rep(survs, N), nrow=length(btime))
   }, list = {
       if (!inherits(pred[[1]], "survfit")) stop("pred is not a list of survfit objects") 
       if (length(pred) != N) stop("pred must be of length(time)")
       pred <- pred[ot]
       survs <-  matrix(unlist(lapply(pred, getsurv, times = btime)),
                               nrow=length(btime), ncol=N)
   }, vector = {
       if (length(pred) != N) stop("pred must be of length(time)")
       if (length(btime) == 1) stop("cannot compute integrated Brier score with pred")
       survs <- pred[ot]
   }, matrix = {
       # <FIXME>
     print(c(length(btime), N))
     print(dim(pred))
       if (all(dim(pred) == c(length(btime), N))){
           survs <- pred[,ot]
       }else{
           stop("wrong dimensions of pred")}
       # </FIXME>
   })
  
   if (is.null(survs)) stop("unknown type of pred")

   # reverse Kaplan-Meier: estimate censoring distribution

   ### deal with ties
   hatcdist <- prodlim(Surv(time, cens) ~ 1,reverse = TRUE)
   # hatcdist <- survfit(Surv(time, 1 - cens) ~ 1)
   # csurv <- getsurv(hatcdist, time)
   # csurv[csurv == 0] <- Inf

   bsc <- rep(0, length(btime))

   # compute Lebesque-integrated Brier score

   if (length(btime) > 1) {
       csurv <- predict(hatcdist, times = btime, type = "surv")
       #csurv <- predict(hatcdist, times = time, type = "surv")
       csurv[csurv == 0] <- Inf

       for (j in 1:length(btime)) {
           help1 <- as.integer(time <= btime[j] & cens == 1)
           help2 <- as.integer(time > btime[j])
           bsc[j] <-  mean(((0 - survs[j])^2) *help1*(1/csurv[j]) +
                           ((1-survs[j])^2)*help2*(1/csurv[j]), na.rm=T)
       }
      # plot(time, bsc,col="green", lty="l")
       ### apply trapezoid rule
       idx <- 2:length(btime)
       RET <- diff(btime) %*% ((bsc[idx - 1] + bsc[idx]) / 2)
       RET <- RET / diff(range(btime))

       ### previously was
       #diffs <- c(btime[1], btime[2:length(btime)] -
       #                     btime[1:(length(btime)-1)])
       #RET <- sum(diffs*bsc)/max(btime)
       names(RET) <- "integrated Brier score"
       attr(RET, "time") <- range(btime)

   # compute Brier score at one single time `btime'
      RET <- list(ibs=RET, bsc=bsc, time=time)
   } else {
       help1 <- as.integer(time <= btime & cens == 1)
       help2 <- as.integer(time > btime)
       cs <- predict(hatcdist, times=btime, type = "surv")
       ### cs <- getsurv(hatcdist, btime)
       if (cs == 0) cs <- Inf
       RET <-  mean((0 - survs)^2*help1*(1/cs) +
                    (1-survs)^2*help2*(1/cs), rm.na=T)
       names(RET) <- "Brier score"
       attr(RET, "time") <- btime
   }
   RET
}



# library(JM)
# ?aucJM

# # install.packages('nftbart')
# library(pec)



# plot(survfit(itx~1))

pred_surv_marg <- function(shape_params, scale_params, mix.prob,data_wide, x) {
  N <- nrow(data_wide)
  k <- nrow(as.matrix(shape_params))

  # f_td <- matrix(0,nrow = N, ncol = k)
  # w <- model.matrix(~ -1+X3+X4+X5, data=data_wide)
  delta <- data_wide$delta
  # initial density function for Weibull
  pred <- matrix(0, nrow = N, ncol = k)
  library(tidyverse)
  dat <- cbind(data_wide,mix.prob)
  colnames(dat) <- c("ID", "event_time", "delta", "V1", "V2")
   dat <- dat %>% arrange((event_time))
   mix.prob <- dat |> select(c(V1,V2)) |> as.matrix()
  # ti <- seq(0, max(data_wide$event_time), length.out = 1000)
  ti <- dat$event_time

  for (j in 1:k) {
    # j <- 1
    shape <- shape_params[j]
    scale <- scale_params[j]
    # gammas <- gammas_matrix[j, ]
    delta <- delta
    # fixed effects in survival component
    # wxi <- as.matrix(w)
    
    pred[,j]<-mix.prob[,j]*pweibull(ti,shape, scale, lower.tail = F)
  }
  pred <- (rowSums(pred))
  pred <- sort(pred, decreasing = T)
  # plot(ti, pred, type = "l")
  # lines(survfit(Surv(ti,delta)~1))
  return(cbind(ti,pred))
}
