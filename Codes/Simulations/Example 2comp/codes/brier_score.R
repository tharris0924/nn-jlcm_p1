get_lcmm_brier_score <- function(model, train_data, test_data, survvars, classmb) {
 #* determine class belongings
  test_data <- data.table(test_data)
  train_data <- data.table(train_data)

 obj <- model

  # mbvars
  # classmb <- paste0("~",mbvars, sep="")
 nclasses <- ncol(obj$pprob) - 2
 coefs <- obj$best
 coefend <- min(which(grepl("Weibull", names(coefs)))) - 1
 coefs_multilogit <- matrix(coefs[1:coefend], nrow = nclasses - 1)
 mbvars <- labels(terms(classmb)) #* get the names of the variables in the model
 # fits <- TRUE
 tmpX <- cbind(1, test_data[, ..mbvars]) #* create a design matrix
 linearval <- as.matrix(tmpX) %*% t(coefs_multilogit) #* get the linear predictions
 test_class <- t(apply(linearval, 1, function(x) {
   exps <- exp(c(x, 0))
   exps / sum(exps)
 })) #* calculate the responsibilities

 #* determine survival curves
 times <- model$predSurv[, 1]
#  ntimes <- length(times)
 ntimes <- length(times)
 # ISE <- 0
 nclasses <- ncol(model$pprob) - 2 # CORRECT
 coefs <- model$best
 coefstart <- max(which(grepl("Weibull", names(coefs)))) + 1
 pred_slopes <- matrix(coefs[coefstart:(coefstart + nclasses - 1)], nrow = nclasses, ncol = length(survvars), byrow = TRUE)
 predclass <- test_class
 if (length(dim(predclass)) == 0) {
   avg <- FALSE
 } else {
   avg <- TRUE
 }
 data <- train_data
 Shat <- NULL
 for (x in c(1:nrow(data))) {
    x <- 1
   if (!avg) {
     tmpc <- predclass[x]
     tmpebx <- exp(sum(pred_slopes[tmpc, ] * data[x, c("X3", "X4", "X5")]))
     Shat <- exp(-tmpebx * model$predSurv[, paste0("event1.CumRiskFct", tmpc)])
   } else {
     Shat_raw <- matrix(0, ncol = nclasses, nrow = ntimes)
          # tmpc <- predclass[x]

     for (tmpc in c(1:nclasses)) {
        tmpc <- 1
      #  tmpebx <- exp(sum(pred_slopes[tmpc, ] * data[x, ..survvars])) 
       Shat_raw[, tmpc] <- exp(-model$predSurv[, paste0("event1.CumRiskFct", tmpc)]) #* multiplying with the baseline cumulative hazard function estimated from the model:
       #* NOTE, this is the same for all i. What makes this change is when it is multiplied with the mixing proportions below
     }
     Shat <- rbind(Shat, c(Shat_raw %*% predclass[x, ]))
   }
 }

 #* calculate integrated Brier score
 evaltimes <- seq(min(data$event_time), max(data$event_time), length.out = 100)
#  torm <- c(which(evaltimes < min(test_data$agedem)), which(evaltimes > max(test_data$agedem)))

 SHAT_k <- Shat
#  SHAT_k <- SHAT_k[, -torm]
 brier_surv <- Surv(test_data$event_time, test_data$delta)
 IBS <- sbrier(brier_surv, t(SHAT_k), evaltimes)
 return(IBS)
}


test_data <- data_wide
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
  k <- nrow(as.matrix(shape_params))

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


get_km_brier_score <- function(test_data) {
 #* determine class belongings
  test_data <- data.table(test_data)
  test_data$rownames <- c(1:nrow(test_data))
 rowtokeep <- test_data[, rownames[1], by = ID][, V1]
   test_data <- test_data[rowtokeep,]
  
    # pred <- matrix(0, nrow = N, ncol = k)
  library(tidyverse)
  dat <- cbind(test_data)
   dat <- dat %>% arrange((event_time))
  #  mix.prob <- dat |> select(c(V1,V2)) |> as.matrix()
  # ti <- seq(0, max(test_data$event_time), length.out = 100)
   ti <- c(dat$event_time)
  # mbvars
  # classmb <- paste0("~",mbvars, sep="")
  
      N <- nrow(test_data)
  k <- nrow(as.matrix(shape_params))

  # initial density function for Weibull
  pred <- matrix(0,nrow=length(ti), ncol=k )
times <- dat$event_time
deltas <- dat$delta
itx2 <- survfit(Surv(times, deltas)~1)$surv
itx <- Surv(times, deltas)
lines(survfit(itx~1), col="red")
      # ti <- test_data$event_time
      # status <- result$data_wide$delta
  

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
dd <- matrix(rep(itx2, 1000), nrow=1000)
 brier_surv <- Surv(test_data$event_time, test_data$delta)
 IBS <- sbrier(obj = brier_surv,pred =  (itx2),btime = evaltimes)
  return(IBS)
}





result$data_long$x <- result$x[result$data_long$ID]
get_lcmm_brier_score(model=jlcmm_fit, train_data = result$data_long, test_data = result$data_long, survvars = NULL, classmb=~x)
 # SHAT[cvid == k, , ] <- SHAT_k
(jlcmm_fit$pprob[,3:4]) <- c("V1","V2")

obj1 <- get_nnem_brier_score(shape_params=shape_est[1,], scale_params=scale_est[1,], mix.prob =jlcmm_fit$pprob[,3:4] , test_data=result$data_wide)
obj2 <- get_nnem_brier_score(shape_params=shape_est[1,], scale_params=scale_est[1,], mix.prob = comp.p, test_data=result$data_wide)
# poly(data$age_init, degree = 2, raw = TRUE)
a <- sbrier(Surv(times, deltas), )

plot(obj1$IBS$time, obj1$IBS$bsc, type="l")
lines(obj2$IBS$time, obj2$IBS$bsc, col="red")
lines(IBS$time, IBS$bsc, col="green")

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



library(JM)
?aucJM

# install.packages('nftbart')
library(pec)



plot(survfit(itx~1))

pred_surv_marg <- function(shape_params, scale_params, mix.prob,data_wide) {
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

est_surv <- pred_surv_marg(shape_params=shape_est[1,], scale_params=scale_est[1,], mix.prob = pij[[1]], data_wide=result$data_wide)
plot(est_surv[,1], est_surv[,2], type="l")
times <- result$data_wide$event_time
deltas <- result$data_wide$delta
itx2 <- survfit(Surv(times, deltas)~1)
itx <- Surv(times, deltas)

 brier_surv <- Surv(times, deltas)
 IBS <- sbrier(obj = brier_surv,pred =  (itx2),btime = times)
lines(survfit(itx~1), col="red")
qt_times <- quantile(times)
bsc <- numeric(length(qt_times))


for(i in 1:5)
    bsc[i] <- sbrier(itx, est_surv[,2], btime=qt_times[i])

data(cancer, package = "survival")
lung <- cancer
lung$status <- lung$status - 1  # Convert to 0/1

# Fit survival models
cox_model <- coxph(Surv(times, deltas) ~ 1, x=TRUE)
km_model <- prodlim(itx~ 1)

library(pec)
# Calculate Brier score at specific time points
brier_scores <- pec(list("NNEM" = cbind((Shat), 0)),
                   formula = Surv(event_time, delta) ~ 1,
                   data = result$data_wide)

print(brier_scores)
pred <- est_surv[,1]
something <- matrix(rep(pred,1000), nrow=1000,ncol=1000)
# Plot prediction error curves
plot(brier_scores)

# Get integrated Brier score
crps(brier_scores)

n <- 100
time <- rpois(n, 20)
cens <- rep(1, n)

# checks, Graf et al. page 2536, no censoring at all!
# no information: \pi(t) = 0.5 

a <- ipred::sbrier(Surv(time, cens), as.matrix(rbind(rep(0.5, n),rep(0.5, n))), time[c(50,51)])
stopifnot(all.equal(cleans(a),0.25))




est_ibs <- sbrier(obj=itx,pred=itx2$surv, btime=1)
obj <- sbrier(obj=itx,pred=itx2$surv)
obj <- sbrier(obj=itx,pred=est_surv[,2])
plot(obj$time, obj$bsc)
Cindex(Surv(times, deltas), est_surv[,2])
nftbart::Cindex(est_surv[,2], times, deltas)

preds <- est_surv
times <- result$data_wide$event_time
deltas <- result$data_wide$delta
concordant_pairs<-0
discordant_pairs<-0

for(i in 1:length(times)){
    for(j in 1:length(times)){
        if(i!=j){
            concordant_pairs <- concordant_pairs+(preds[i]>preds[j])*(times[i]>times[j])*deltas[j]
            discordant_pairs <- discordant_pairs+(preds[i]<preds[j])*(times[i]>times[j])*deltas[j]

        }
    }
}

print(concordant_pairs)


c_index <- concordant_pairs/(concordant_pairs+discordant_pairs)
print(c_index)
library("survival")
library(ipred)
data("DLBCL", package = "ipred")
smod <- Surv(DLBCL$time, DLBCL$cens)

KM <- survfit(smod ~ 1)
# integrated Brier score up to max(DLBCL$time)
sbrier(smod, KM)
range(smod[,1])

library(yardstick)

# Install if needed: install.packages("pec")
library(pec)
library(survival)

# Example with lung cancer data
data(cancer)
cancer$status <- cancer$status - 1  # Convert to 0/1 (0=censored, 1=event)

# Fit survival models
fit_cox <- (coxph(Surv(time, status) ~ age + sex + ph.ecog, data = cancer, x=T))
fit_km <- survfit(Surv(time, status) ~ 1, data = cancer)
plot(summary(fit_cox)$surv)
time_points <- seq(100, 1000, 100)
  for (i in seq_along(time_points)) {
    t <- time_points[1]
    
    # Get survival probabilities at time t
    surv_pred <- summary(survfit(fit_cox, newdata = cancer), times = t)$surv
    
    # Calculate Brier score at time t
    brier_scores[i] <- brier_survival(surv_pred, t, cancer$time, cancer$status)
  }

# Calculate Brier score at specific time points
brier_result <- pec(list("Cox" = fit_cox),
                    formula = Surv(time, status) ~ 1,
                    data = cancer,
                    times = c(100, 200, 365, 500))

print(brier_result)
plot(brier_result)


predSurv <- function(shape_params, scale_params, mix.prob,data_wide) {
  N <- nrow(data_wide)
  k <- nrow(as.matrix(shape_params))

  # f_td <- matrix(0,nrow = N, ncol = k)
  # w <- model.matrix(~ -1+X3+X4+X5, data=data_wide)
  delta <- data_wide$delta
  # initial density function for Weibull
  pred <- matrix(0,nrow=100, ncol=k )
  # library(tidyverse)
  # dat <- cbind(data_wide,mix.prob)
  # dat <- dat %>% arrange((event_time))
  # mix.prob <- dat |> select(c(V1,V2))
      ti <- seq(0, max(result$data_wide$event_time), length.out = 100)
      status <- result$data_wide$delta
  for (j in 1:k) {
    # j <- 1
    shape <- shape_params[j]
    scale <- scale_params[j]
    # gammas <- gammas_matrix[j, ]
    delta <- delta
    # fixed effects in survival component
    # wxi <- as.matrix(w)
    pred[,j] <- pweibull(ti,shape, scale, lower.tail = F)
  }
  plot(ti, pred[,1], type="l")
  lines(ti, pred[,2])
  # lines(survfit(Surv(ti,status)~1, weights=mix.prob[,2]), bty="l")

  return(pred)
}

# nftbart::Cindex()
model=jlcmm_fit
predclass=pij[[1]]

data=result$data_wide
x <- result$x
ISE1 <- get_lcmm_ISE(model=jlcmm_fit, data=result$data_wide,groups=result$x)
get_lcmm_ISE <- function(model, data,groups, true_shape, true_scale, true_pi){
    times <- model$predSurv[,1]
    ntimes <- length(times)
    ISE <- 0

  
    nclasses <- ncol(model$pprob)-2
    coefs <- model$best
    coefstart <- max(which(grepl('Weibull',names(coefs))))+1
   coefend <- min(which(grepl("Weibull", names(coefs)))) - 1

    # pred_slopes <- matrix(coefs[coefstart:(coefstart+nclasses*3-1)],nrow=nclasses,ncol=3)
     coefs_multilogit <- matrix(coefs[1:coefend], nrow = nclasses - 1)
#  mbvars <- labels(terms(classmb)) #* get the names of the variables in the model
 # fits <- TRUE
 tmpX <- cbind(1, groups) #* create a design matrix
 linearval <- as.matrix(tmpX) %*% t(coefs_multilogit) #* get the linear predictions
 test_class <- t(apply(linearval, 1, function(x) {
   exps <- exp(c(x, 0))
   exps / sum(exps)
 })) #* calculate the responsibilities
 predclass <- test_class
    if (length(dim(predclass))==0){ avg <- FALSE} else {avg <- TRUE}
    for (x in c(1:nrow(data))){
        if (!avg){
            tmpc <- predclass[x]
            tmpebx <- exp(sum(pred_slopes[tmpc,] * data[x,c('X3','X4','X5')]))
            Shat <- exp(-tmpebx*mod$predSurv[,paste0('event1.CumRiskFct',tmpc)])
        } else {
            Shat_raw <- matrix(0,ncol=nclasses,nrow=ntimes)
            for (tmpc in c(1:nclasses)){
                # tmpebx <- exp(sum(pred_slopes[tmpc,] * data[x,c('X3','X4','X5')]))
                 Shat_raw[, tmpc] <- exp(-model$predSurv[, paste0("event1.CumRiskFct", tmpc)])
            }
            Shat<- c(Shat_raw %*% predclass[x,]) #* individual survival curve at time 
          # Shat <- rbind(Shat, c(Shat_raw %*% predclass[x, ]))
        }
        # Shat <- pred_surv_marg(shape_params = result$theta$shape, result$theta$scale, t(as.matrix(result$theta$pj[x,])), result$data_wide)
        Strue <- pred_surv_marg(shape_params = true_shape, true_scale, t(as.matrix(true_pi)), result$data_wide)


        scores <- (Shat - Strue)^2
        SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(times)) / diff(range(times))
        ISE <- ISE + SE
    }
    ISE <- ISE / nrow(result$data_wide)
    return (ISE)
}
mix.prob <- pij[[1]]
shape_est <- 
get_nnem_ISE <- function(result, shape_est, scale_est, mix.prob,theta_orig, true_pi){
  ISE <- 0
      times <- seq(0,max(result$data_wide$event_time), length.out = 100)
    ntimes <- length(times)
      for (x in c(1:nrow(result$data_wide))){
        Shat <- pred_surv_marg(shape_params = shape_est, scale_est, t(as.matrix(mix.prob[x,])), result$data_wide)
        Strue <- pred_surv_marg(shape_params = theta_orig$shape, scale_params=theta_orig$scale,mix.prob=t(as.matrix(true_pi[x,])), data_wide=result$data_wide)
        scores <- (Shat[,2] - Strue[,2])^2
        SE <- sum(0.5*(scores[-1]+scores[-ntimes]) * diff(times)) / diff(range(times))
        ISE <- ISE + SE
    }
    ISE <- ISE / nrow(result$data_wide)
  return(ISE)
}


