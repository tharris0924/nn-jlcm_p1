k <- 2
library(lcmm)
library(tictoc)
# init <- Jointlcmm(fixed = y ~ measure_time,
#   #  mixture=~1,
#   # classmb = ~ X1 + X2,
#   random = ~1,
#    subject='ID',
#    survival = Surv(event_time, delta) ~ 1,
#    hazard="Weibull",hazardtype="Specific", ng=1, data=data, logscale = T)
#  summary(init)
data$x <- x[data$ID]
# n_cores <- parallel::detectCores() - 1
init <- matrix(nrow = 12)
init[1:2] <- c(0,0) 
init[3:6] <- exp(init[3:6])
scales_est <- init[c(5,3)]^(-1/init[c(6,4)])
init[c(5,3)] <- log(scales^(-shapes))
init[c(6,4)] <- log(shapes)
init[c(8,7)]<-betas_matrix[,1]
init[c(10,9)] <- betas_matrix[,2]
init[c(11)] <- D
init[c(12)] <- sd_e
# init <- as.vector(init)
jlcmm_fit <- tryCatch({Jointlcmm(
  fixed = y ~ measure_time,
  mixture=~measure_time,
  classmb = ~ x,
  random = ~1,
  subject='ID',
  survival = Surv(event_time, delta) ~ 1,
  hazard="Weibull",
  hazardtype="Specific", 
  ng=2, 
  data=data, nwg=F,logscale = T, B=as.vector(init)
)}, error=function(e){NULL})
summary(jlcmm_fit)
# fixef(jlcmm_fit)
init <- jlcmm_fit$best
