# k <- 2
library(lcmm)
library(tictoc)
init <- Jointlcmm(fixed = y ~ measure_time + measure_time2,
  #  mixture=~1,
  # classmb = ~ X1 + X2,
  random = ~measure_time + measure_time2,
   subject='ID',
   survival = Surv(event_time, delta) ~ 1,
   hazard="Weibull",hazardtype="Specific", ng=1, data=data, logscale = T)
 summary(init)
data$x <- x[data$ID]
# # n_cores <- parallel::detectCores() - 1
# # init[1:2] <- c(0,0) 
# # init[3:6] <- exp(init[3:6])
# # scales_est <- init[c(5,3)]^(-1/init[c(6,4)])
# init <- estimates(mj4b)
# scales_est <- init[c(10,8,6,4)]^(-1/init[c(11,9,7,5)])
# shapes_est <- init[c(11,9,7,5)]
# betas_matrix <- round(matrix(init[14:25], ncol = 3))
# init <- matrix(0,nrow=(19+8+3))
# init[1:3] <- 0
# init[c(10,8,6,4)] <- log(scales^(-shapes))
# init[c(11,9,7,5)] <- log(shapes)
# init[c(15,14,13,12)]<-betas_matrix[,1]
# init[c(19,18,17,16)] <- betas_matrix[,2]
# init[c(23,22,21,20)] <- betas_matrix[,3]

# init[c(24:29)] <- vech((chol(D_g_og)))
# init[c(30)] <- 10
# init<-c(0,0,0,init)
# # init[22] <- 1
# jlcmm_fit <- tryCatch({Jointlcmm(
#   fixed = y ~ measure_time + measure_time2,
#   mixture=~measure_time + measure_time2,
#   classmb = ~ x,
#   random = ~measure_time + measure_time2,
#   subject='ID',
#   survival = Surv(event_time, delta) ~ 1,
#   hazard="Weibull",
#   hazardtype="Specific", 
#   ng=4, 
#   data=data, nwg=F,logscale = T, B=random(init), nproc=detectCores()-1
# )}, error=function(e){NULL})
# summary(jlcmm_fit)
# fixef(jlcmm_fit)
jlcmm_fit<-tryCatch({gridsearch(rep = 15, maxiter = 10, minit = init, Jointlcmm( y ~ measure_time + measure_time2,
   mixture =~ measure_time + measure_time2, random =~ measure_time + measure_time2, 
   survival = Surv(event_time, delta) ~ 1,
    hazard = "Weibull", subject = 'ID', data = data, ng = 4,verbose=F, logscale = T),cl=detectCores()-1)},error=function(e){NULL})
summary(jlcmm_fit)
# init <- jlcmm_fit$best
# init[1:4] <- 0
# init[5:10] <- exp(init[5:10])
# scales_est <- init[c(9,7,5)]^(-1/init[c(10,8,6)])

