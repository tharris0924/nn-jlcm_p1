# k <- 2
library(lcmm)
library(tictoc)
init <- Jointlcmm(fixed = y ~ measure_time,
  #  mixture=~1,
  # classmb = ~ X1 + X2,
  random = ~measure_time,
   subject='ID',
   survival = Surv(event_time, delta) ~ 1,
   hazard="Weibull",hazardtype="Specific", ng=1, data=data, logscale = T)
 summary(init)
data$x <- x[data$ID]
# n_cores <- parallel::detectCores() - 1
# init[1:2] <- c(0,0) 
# init[3:6] <- exp(init[3:6])
# scales_est <- init[c(5,3)]^(-1/init[c(6,4)])
init <- matrix(0,nrow=22)
init[1:4] <- 0
init[c(9,7,5)] <- log(scales^(-shapes))
init[c(10,8,6)] <- log(shapes)
init[c(13,12,11)]<-betas_matrix[,1]
init[c(16,15,14)] <- betas_matrix[,2]
init[c(17:19)] <- vech(D)
init[c(20:21)] <- omegas
init[22] <- 1
jlcmm_fit <- tryCatch({Jointlcmm(
  fixed = y ~ measure_time,
  mixture=~measure_time,
  classmb = ~ x,
  random = ~measure_time,
  subject='ID',
  survival = Surv(event_time, delta) ~ 1,
  hazard="Weibull",
  hazardtype="Specific", 
  ng=3, 
  data=data, nwg=T,logscale = T, B=as.vector(init)
)}, error=function(e){NULL})
summary(jlcmm_fit)
# fixef(jlcmm_fit)

# init <- jlcmm_fit$best
# init[1:4] <- 0
# init[5:10] <- exp(init[5:10])
# scales_est <- init[c(9,7,5)]^(-1/init[c(10,8,6)])

