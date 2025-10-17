
library("lcmm")
library("NormPsy")
library('data.table')
library(tictoc)
library('survival')
paquidS <- paquid[which(paquid$agedem> paquid$age_init),]
paquidS$normMMSE <-normMMSE(paquidS$MMSE)
paquidS$age65 <- (paquidS$age - 65)/10
seed <- 1; n <- 1000; censor_rate <- 0.05


#! test - train split
library(jlctree)
library(partykit)
library(tictoc)
# source("utils.R")
library("NormPsy")
library(lcmm)
library(doParallel)
library(data.table)

# Construct data
paquid$normMMSE <- normMMSE(paquid$MMSE)
paquid$age65 <- (paquid$age - 65) / 10
head(paquid)
paquidS <- paquid[paquid$agedem > paquid$age_init & paquid$age <= paquid$agedem, ]

paquidS2 <- data.table(paquidS)
paquidS2$age <- paquidS2[,
  {
    if (.N == 1) {
      age_init
    } else {
      c(age_init[1], age[c(1:(.N - 1))])
    }
  },
  by = ID
][, V1]

temp2 <- subset(paquidS2, select = c(ID, age_init, agedem, dem, male, CEP)) # baseline
data_wide <- unique(temp2)

temp <- subset(paquidS2, select = c(ID, age_init, agedem, dem, male)) # baseline
temp <- unique(temp)
data <- tmerge(temp, temp, id = ID, tstart = age_init, tstop = agedem, death = event(agedem, dem)) # set range
data <- tmerge(data, paquidS2,
  id = ID,
  age65 = tdc(age, age65), CEP = tdc(age, CEP),
  normMMSE = tdc(age, normMMSE),
  BVRT = tdc(age, BVRT),
  IST = tdc(age, IST),
  HIER = tdc(age, HIER),
  CESD = tdc(age, CESD)
)

data <- data.table(data)
data <- data[!is.na(normMMSE) & !is.na(BVRT) & !is.na(IST) & !is.na(HIER) & !is.na(CESD)]

data_tmp <- data.table(data)
data_tmp <- data_tmp[, list(
  age65_one = rep(age65[1], .N),
  BVRT_one = rep(BVRT[1], .N),
  IST_one = rep(IST[1], .N),
  HIER_one = rep(HIER[1], .N),
  CESD_one = rep(CESD[1], .N)
), by = ID]
data$age65_one <- data_tmp$age65_one
data$BVRT_one <- data_tmp$BVRT_one
data$IST_one <- data_tmp$IST_one
data$HIER_one <- data_tmp$HIER_one
data$CESD_one <- data_tmp$CESD_one


data$age65_sq <- data$age65^2
data <- data.table(data)
ndata <- nrow(data)
nfolds <- 5
n_per_fold <- ceiling(ndata / nfolds) * c(1:nfolds)
N_by_id <- data[, .N, by = ID]
set.seed(0)
id_shuffle <- sample(seq_len(nrow(N_by_id)), replace = FALSE)
N_by_id <- N_by_id[id_shuffle, ]
N_by_id$cumN <- cumsum(N_by_id$N)
N_by_id$cvid <- apply(array(N_by_id$cumN), 1, function(x) {
  sum(x > n_per_fold) + 1
})
data <- merge(data, N_by_id[, list(ID, cvid)], by = "ID")
cvid <- data$cvid

temp2 <- subset(data, select = c(ID, agedem, dem, male, CEP, cvid)) # baseline
data_wide <- unique(temp2)
data_long <- data

cv <- 2
train_data_long <- data[cvid != cv, ]
test_data_long <- data[cvid == cv, ]
train_data_wide <- data_wide[cvid != cv, ]
test_data_wide <- data_wide[cvid == cv, ]

# library(lcmm)
# data(paquid)
m <- 2

?ggsurvfit::ggsurvfit
m <- 2
# computation of the normalized MMSE
# paquid$MMSEnorm <- normMMSE(paquid$MMSE)

# histogram of these data
# par(mfrow=c(1,2))
# hist(paquid$MMSE,breaks=seq(-0.5,30.5,1),col=2,main="crude MMSE")
# hist(paquid$MMSEnorm,col=3,main="normalized MMSE")
# paquidS$logCESD <- (log(paquidS$CESD))
mj1 <- jlcmm(fixed = normMMSE ~ age65 +I(age65^2) + CEP,
             random =~ age65+I(age65^2),
             survival = Surv(age_init, agedem,dem)~CEP+male, hazard = "Weibull",
             subject = 'ID', data = train_data_long, verbose = F, logscale = T, maxiter=250)
# summary(survival::Surv(paquidS$agedem, paquidS$dem))
summary(mj1)  

# mj1$best[1:2]^2

# mj2b <- jlcmm(fixed =normMMSE ~ age65 + I(age65^2) + CEP,
#                         classmb = ~ CEP+male,
#                         random = ~ age65 + I(age65^2), 
#                         mixture = ~ age65 + I(age65^2), 
#                         ng = 2,
#                         survival = Surv(age_init, agedem, dem) ~ CEP + male, 
#                         hazard = "Weibull",
#                         subject = 'ID', data = train_data_long,
#                         B = random(mj1), logscale = T, nproc=19, maxiter=250)
# summary(mj2b)
# save(mj2b,file= "NNEM/Application/init2model_application.RData")
# # sqrt(diag(vcov(mj2b)))

# mj3b <- jlcmm(fixed =normMMSE ~ age65 + I(age65^2) + CEP,
#                         classmb = ~ CEP+male,
#                         random = ~ age65 + I(age65^2), 
#                         mixture = ~ age65 + I(age65^2), 
#                         ng = 3,
#                         survival = Surv(age_init, agedem, dem) ~ CEP + male, 
#                         hazard = "Weibull",
#                         subject = 'ID', data = train_data_long,
#                         B = random(mj1), logscale = T, nproc=19, maxiter = 250)

# # # tic()
# summary(mj3b)
# save(mj3b,file= "NNEM/Application/init3model_application.RData")

# mj4b <-jlcmm(fixed =normMMSE ~ age65 + I(age65^2) + CEP,
#                         classmb = ~ CEP+male,
#                         random = ~ age65 + I(age65^2), 
#                         mixture = ~ age65 + I(age65^2), 
#                         ng = 4,
#                         survival = Surv(age_init, agedem, dem) ~ CEP + male, 
#                         hazard = "Weibull",
#                         subject = 'ID', data = train_data_long,
#                         B = random(mj1), logscale = T, nproc=19, maxiter = 250)
                  
# summary(mj4b)
# save(mj4b,file= "NNEM/Application/init4model_application.RData")


# summarytable(mj2b, mj3b, mj4b, which=c("loglik","AIC", "BIC", "scoretest", "entropy"))
# # ST <- c(mj1$scoretest[1],mj2b$scoretest[1],mj3b$scoretest[1],mj4b$scoretest[1]) 
# # plot(ST ~ c(1:4),type="l",bty='n',main="Score test", xlab="number of classes", ylab="")
# # ??jlcmm
# # summary(mj4b)
# # plot(mj4b)
# # toc()

# # library(joineRML)
# # ?joineRML::mjoint()
# library(JM)

#  lme_model_fit <- lme(normMMSE ~ age65+ I(age65^2)+CEP, random = ~ age65 + I(age65^2)|ID, data = train_data_long)
# surv_model_fit <- coxph(Surv(agedem, dem)~CEP+male, data=train_data_wide, x=TRUE, model=TRUE)
# fit1 <- jointModel(lmeObject = lme_model_fit, surv_model_fit, timeVar = "age65")
# summary(fit1)
# # ?jointModel

# preds <- survfitJM(fit1, train_data_long, idVar = "ID", simulate=F, survTimes = sort(train_data_wide$agedem))

# # save(preds, file="predSurvJM.RData")

# # time_points <- preds$summaries[[1]][,1]

# # # Combine all predSurv columns column-wise
# # surv_matrix <- sapply(preds$summaries, function(x) {
# #   # handle if x is a list, data.frame, or atomic
# #   if (is.data.frame(x)) {
# #     x$predSurv
# #   } else if (is.list(x) && !is.null(x$predSurv)) {
# #     x$predSurv
# #   } else if (is.matrix(x) && "predSurv" %in% colnames(x)) {
# #     x[, "predSurv"]
# #   } else {
# #     stop("Each element must contain a 'predSurv' component.")
# #   }
# # })
library('ggsurvfit')
survfit2(Surv(agedem, dem) ~ 1, data = train_data_wide) %>%
  ggsurvfit() +
  add_censor_mark() +
  add_confidence_interval() +
  add_quantile() +
  add_risktable() +
  scale_ggsurvfit(x_scales = list(limits=c(65,100)))+
  theme_ggsurvfit_KMunicate()

# # Transpose so that rows = subjects, columns = time points
# surv_matrix <- t(surv_matrix)

# # Assign row and column names
# rownames(surv_matrix) <- names(preds$summaries)  # subject IDs
# colnames(surv_matrix) <- round(time_points, 2)
# brier_surv <- ret$brier_surv
# it<- IBS(object=brier_surv, sp_matrix=surv_matrix, IBSrange=sort(data_wide$event_time), n_cores=19,plot=TRUE, return_brier = TRUE)
# it_srem <- it
# save(it_srem, file="IBS_srem.RData")


#  lme_model_fit <- lme(normMMSE ~ age65+ I(age65^2)+CEP, random = ~ age65 + I(age65^2)|ID, data = train_data_long)
# surv_model_fit <- coxph(Surv(agedem, dem)~CEP+male, data=train_data_wide, x=TRUE, model=TRUE)
# fit1 <- jointModel(lmeObject = lme_model_fit, surv_model_fit, timeVar = "age65")
# summary(fit1)
# ?jointModel

# preds_test <- survfitJM(fit1, test_data_long, idVar = "ID", simulate=F, survTimes = sort(test_data_wide$event_time))
# save(preds_test, file="test_predSurvJM.RData")


# time_points <- preds_test$summaries[[1]][,1]

# # Combine all predSurv columns column-wise
# surv_matrix_test <- sapply(preds_test$summaries, function(x) {
#   # handle if x is a list, data.frame, or atomic
#   if (is.data.frame(x)) {
#     x$predSurv
#   } else if (is.list(x) && !is.null(x$predSurv)) {
#     x$predSurv
#   } else if (is.matrix(x) && "predSurv" %in% colnames(x)) {
#     x[, "predSurv"]
#   } else {
#     stop("Each element must contain a 'predSurv' component.")
#   }
# })


# # Transpose so that rows = subjects, columns = time points
# surv_matrix_test <- t(surv_matrix_test)

# # Assign row and column names
# rownames(surv_matrix_test) <- names(preds_test$summaries)  # subject IDs
# colnames(surv_matrix_test) <- round(time_points, 2)
# brier_surv <- Surv(test_data_wide$event_time, test_data_wide$delta)
# # brier_surv <- ret$brier_surv
# it<- IBS(object=brier_surv, sp_matrix=surv_matrix_test, IBSrange=sort(test_data_wide$event_time), n_cores=19,plot=TRUE, return_brier = TRUE)
# it_srem_test <- it
# save(it_srem_test, file="IBS_srem_test.RData")

# # Inspect
# surv_matrix[1:3, 1:5]

# preds$summaries
# fit1 <- mjoint(formLongFixed = normMMSE ~ age65 + I(age65^2) + CEP,
#     formLongRandom = ~ age65 + I(age65^2)|ID,
#     formSurv = Surv(agedem, dem) ~ CEP + male,
#     survData = train_data_wide,
#     data = train_data_long,
#     timeVar = "age65") # controls for illustration only
# summary(fit1)
# plot(mj3b,which="fit",var.time="age65" ,bty="l",xlab = "Decades from 65", ylab = "normMMSE", marg=FALSE)
# # load("NNEM/Application/init2model_application.RData")
# # summary(mj4b)
# # library(lme4)
# # lmer(normMMSE ~ age65 + I(age65^2) + CEP + (age65|ID),paquidS)
 init <- estimates(mj1, cholesky=F)
# D <- VarCovRE(mj4b)

vcovmj1 <- vcov(mj1)
#  matrix(D[,1], nrow=3, ncol=3, byrow = F)


sample_init <- mvnfast::rmvn(m, init, sigma = vcovmj1)
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
# scales_est <- init[c(10,8,6,4)] 

# scales_est <- (init[,1]^2)^(-1/(init[,2]^2))
# shapes_est <- exp(sample_init[,2])

gammas_matrix <- sample_init[1,(max(idWeib_shape)+1):(max(idWeib_shape)+2)]

# betas_matrix <- sample_init[,2:8]

idxs_slope <- which(grepl(c('intercept', "age65", "measure_time2")[2], colnames(sample_init)))

          b1 <- sample_init[,(idxs_slope[1]-1)]
          b2 <- sample_init[,idxs_slope[1]]
          b3 <- sample_init[(idxs_slope[2])]

          df <- df2 <- data.frame(int=b1,mt=b2,mt2=b3)
          # df[2,] <- df2[1,]
          # df[1,] <- df2[3,]
          # df[3,] <- df2[4,]
          # df[4,] <- df[2,]
betas_matrix <- as.matrix(df)
beta <- sample_init[1,(max(idxs_slope)+1)]
# cholesky <- metafor::vec2mat(sample_init[1,9:14], diag=T)
# cholesky[3,1] <-cholesky[2,2] 
# cholesky[2,2] <- cholesky[1,3]
# cholesky[1,3] <- cholesky[3,1]

residual_err <-sample_init[1,"stderr"]

par(mfrow=c(1,1))



datalong <- train_data_long
# data <- datalong[is.finite(rowSums(datalong, na.rm = TRUE)), ]
data <- datalong

nis<-table(datalong$ID)
datalong$delta <- datalong$dem
test_data_long$delta <- test_data_long$dem
datalong$event_time <- datalong$agedem
test_data_long$event_time <- test_data_long$agedem

sim <- data.table(datalong)
datalong$y <- datalong$normMMSE
temp2 <- subset(sim, select = c(ID, event_time, delta, CEP, male)) # baseline
data_wide <- unique(temp2)

sim <- data.table(test_data_long)
test_data_long$y <- test_data_long$normMMSE
temp2 <- subset(sim, select = c(ID, event_time, delta, CEP, male)) # baseline
test_data_wide <- unique(temp2)

n <- nrow(data_wide)

datalong$measure_time2 <- datalong$age65^2
datalong$measure_time <- datalong$age65

# bi <- list()

test_data_long$measure_time2 <- test_data_long$age65^2
test_data_long$measure_time <- test_data_long$age65

# Load required libraries
library(dplyr)
library(ggplot2)


# Set parameters for simulation
# set.seed(123)  # for reproducibility
# seed <- 1;n <- 500; censor_rate <- 0.05
# set.seed(seed)
# # source("NNEM/helper_functions/util.R")

# # example 3: Xue & Yao 2022
# n_samples <- n
# x <- runif(n_samples, -5, 12)
# # X <- cbind(1, x) 
# x <- sort(x)
# true_mixing_prop <- function(x) {
#     pi1 <- (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4)) * as.numeric(x < 3)
#     pi2 <- (1 - (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4))) * as.numeric(x < 0)
#     pi3 <- (2 * exp(-0.1 * (x - 6)^4)) / (1 + exp(-0.1 * (x - 6)^4)) * as.numeric(x >= 3)
#     pi4 <- 1 - (pi1 + pi2+ pi3)
#     return(cbind(pi1, pi2, pi3, pi4))
#   }

# library(dplyr)
pi_init <- matrix(rep(1/m), nrow=nrow(data_wide), ncol=m)

# df <- data.frame(x,pi_init)
# df<-df|>arrange(x)
# plot(df$x,df$pi1, type="l", col="blue")
# lines(df$x, df$pi2)
# lines(df$x,df$pi3, col="green")
# lines(df$x,df$pi4, col="red")
component_assignments <- apply(pi_init, 1, function(p) sample(1:m, 1, prob = p))

time_Y <- paquidS$agedem
delta <- paquidS$dem
hist(time_Y)
# Model parameters
n_subjects <- length(unique(datalong$ID))     # number of individuals
# n_timepoints <-     # number of measure_time points per individual

hist(table(paquidS$ID))
median(table(paquidS$ID))

# 
datalong$y <- datalong$normMMSE

    # # fixef <- c(0.5,2,1)
    # betas_matrix <- matrix(c(0.5,2, 
    #                         3, 2,
    #                         2, 1), ncol = 2, byrow=T)
    # # betas_matrix <- matrix(c(fixef,fixef), nrow = 3, ncol=2)
  
    # # ranefsx <- bi[[pseudo_g]]
    # ranefs <- ranef[datalong$ID,]
    
    # Zb <- ranefs[,1] + ranefs[,2]*datalong$measure_time
    # y <- betas_matrix[,1][pseudo_g]+betas_matrix[,2][pseudo_g]*datalong$measure_time+ Zb + rnorm(nrow(datalong), sd=sd_e) 
library(MASS)
library(matrixcalc)
n_timepoints <- table(data$ID)
# Create design matrix
# Z <- model.matrix(~measure_time, data = data)
# X <- Z
sigma_e <- residual_err
D <- cholesky
  # V_blocks <- list()
    # omega <- props

# mu <- c()
#   for (i in 1:n) {
#     # i <- 1
#     idx <- which(data$ID == i)
#     g <- component_assignments[i]
#     Z_i <- Z[idx, , drop = FALSE]
#     V_i <- Z_i %*% (D) %*% t(Z_i) + (sigma_e^2) * diag(length(idx))
#     X_i <- X[idx, ,drop=FALSE]
#     XBg <- X_i %*% betas_matrix[g, ]
#     V_blocks[[i]] <- V_i
#     mu <- c(mu, XBg)
#   }
# V <- Matrix::bdiag(V_blocks)
# library(mvnfast)
# y <- t(rmvn(1, mu, V))
# datalong$y <- y
# datalong <- paquidS
# datalong$event_time <- datalong$agedem
pseudo_g <- component_assignments[datalong$ID]
# head(datalong,10)
# Clean up intermediate columns (optional): don't need to do this if we generate from a poison
# data <- data %>%
#   select(ID, measure_time, y) %>%
#   arrange(ID, measure_time)

# Display first few rows
# head(datalong, 15)
datalong$g <-pseudo_g 
datalong$measure_time <- datalong$age65
# Summary statistics
cat("Summary of generated data:\n")
cat("Number of subjects:", n_subjects, "\n")
# cat("Number of observations per subject:", mean(as.numeric(n_timepoints)), "\n")
cat("Total observations:", nrow(datalong), "\n")
cat("Mean response:", round(mean(datalong$y, na.rm=T), 2), "\n")
cat("SD of response:", round(sd(datalong$y, na.rm=T), 2), "\n")

# Plot a sample of trajectories
library(dplyr)
sample_subjects <- sample(1:n_subjects, min(50, n_subjects))
data_sample <- datalong %>% filter(ID %in% sample_subjects)

p <- ggplot(data_sample, aes(x = measure_time, y = y, group = ID, color = factor(dem))) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "Time",
       y = "Response (Y)",
       color = "Dementia Diagnosis") +
  theme_minimal()+
  theme(legend.position = "bottom")
print(p)
# Verify the model by fitting it back to the generated data
# library(lme4)
# model_fit <- lmer(y ~ CEP+(measure_time+measure_time2| ID), data = datalong)
# vcov.l,er(model_fit)
# print(summary(model_fit))


  data <- datalong

    m <- length(scale_est)
    # D_g <- D_g_og
    # alphas <- props
    # alphas
    # gammas_matrix <- slopes[1:m,]
    # mu_g <- mu_g



library('data.table')
    beta_o <- beta 

#? compute class specific random effects 
m <- length(shape_est)
bi <- list()
Xg <- model.matrix(~1+ measure_time +(measure_time2), data=datalong)
Xo <- model.matrix(~-1+CEP, datalong)
Zg <-  model.matrix(~1+ measure_time +(measure_time2), data=datalong)
q <- ncol(Zg)
for(j in 1:m){
# j <- 1
    sigma_g <- residual_err # error variances 

    D_g <- D
     # random effect variance-covariance matrices
    beta_g <- betas_matrix[j,] # fixed effects
    beta_o <- beta 

  # i<-1
  bigs<-matrix(nrow = n, ncol = q)
  seq <- data_wide$ID

    for(i in seq_len(length(seq))){

      # i<-2
      # n_i <- nis[i]
      # k <- 1
      ids <- seq[i]
      sim <- data
      subset<-sim[sim$ID==ids,]
      Z_ig <- Zg[sim$ID==ids, , drop = FALSE]
      X_ig <- Xg[sim$ID==ids, ,drop=FALSE]
      X_io <- Xo[sim$ID==ids, ,drop=FALSE]

      p<-ncol(X_ig)
      Xi_beta <- X_ig %*% beta_g + X_io %*% beta_o
      # xbs[i]<-list(Xi_beta)
      q<-ncol(Z_ig)
      Y_i <- subset$y
      # ys<-list(Y_i)
      sigma<-sigma_g
      n_i <- length(subset$ID)
      # 3by2 %*% 
      Vi <- as.matrix(Z_ig%*% D_g%*% t(Z_ig) + sigma_g^2 * diag(n_i))
      big <- D_g %*% t(Z_ig) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
      bigs[i,] <- big
      }
  bi[j] <- list(bigs)
}



library(nnet)
cholesky <- t(chol(D))
tcrossprod(cholesky)
# omegas <- alphas[1:2]
theta_orig <- theta <- list(pj=pi_init, scale=scale_est, shape=shape_est, 
   betas=betas_matrix, beta_o=beta_o, ranef=bi, D=D, cholesky=cholesky, sigma_e=residual_err, gamma_matrix=gammas_matrix)
save(theta_orig, file="theta_orig_k4.RData")



# # Plot true class-specific survival probabilities
# library(ggplot2)
# library(dplyr)

# # Use your existing parameters
# shapes <- round(shapes_est)
# scales <- round(scales_est)

# print("True Weibull parameters:")
# print(data.frame(Class = 1:4, Shape = shapes, Scale = scales))

# # Create time sequence for plotting
# time_seq <- seq(0, 200, by = 1)

# # Function to calculate Weibull survival probability
# weibull_survival <- function(t, shape, scale) {
#   exp(-(t/scale)^shape)
# }

# # Calculate survival probabilities for each class
# survival_data <- data.frame()

# for(class in 1:4) {
#   class_survival <- data.frame(
#     time = time_seq,
#     survival_prob = weibull_survival(time_seq, shapes[class], scales[class]),
#     class = paste("Class", class)
#   )
#   survival_data <- rbind(survival_data, class_survival)
# }

# # Create the plot
# p1 <- ggplot(survival_data, aes(x = time, y = survival_prob, color = class)) +
#   geom_line(size = 1.2) +
#   labs(
#     title = "True Class-Specific Survival Probabilities",
#     subtitle = paste("Weibull distributions with shapes:", paste(shapes, collapse = ", "), 
#                      "and scales:", paste(scales, collapse = ", ")),
#     x = "Time",
#     y = "Survival Probability",
#     color = "Latent Class"
#   ) +
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
#   scale_color_manual(values = c("Class 1" = "#1f77b4", "Class 2" = "#ff7f0e", 
#                                "Class 3" = "#2ca02c", "Class 4" = "#d62728")) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 14, face = "bold"),
#     plot.subtitle = element_text(size = 12),
#     legend.position = "right",
#     panel.grid.minor = element_blank()
#   )

# print(p1)

# # Alternative: Base R plot
# par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
# plot(time_seq, weibull_survival(time_seq, shapes[1], scales[1]), 
#      type = "l", col = "blue", lwd = 2,
#      ylim = c(0, 1), xlim=c(65, 120),xlab = "Time", ylab = "Survival Probability",
#      main = "True Class-Specific Survival Probabilities")

# for(class in 2:4) {
#   lines(time_seq, weibull_survival(time_seq, shapes[class], scales[class]), 
#         col = c("orange", "green", "red")[class-1], lwd = 2)
# }

# legend("topright", 
#        legend = paste("Class", 1:4, "(λ=", scales, ", κ=", shapes, ")"),
#        col = c("blue", "orange", "green", "red"),
#        lwd = 2, cex = 0.8)

# # Additional analysis: median survival times
# median_survival <- numeric(4)
# for(class in 1:4) {
#   # Median survival time for Weibull: scale * (ln(2))^(1/shape)
#   median_survival[class] <- scales[class] * (log(2))^(1/shapes[class])
# }

# print("Median survival times by class:")
# print(data.frame(Class = 1:4, 
#                  Shape = shapes, 
#                  Scale = scales, 
#                  Median_Survival = round(median_survival, 2)))

# # Hazard ratios (relative to Class 1)
# cat("\nHazard characteristics:\n")
# for(class in 1:4) {
#   cat("Class", class, "- Shape:", shapes[class], 
#       "(", ifelse(shapes[class] > 1, "increasing", 
#                   ifelse(shapes[class] == 1, "constant", "decreasing")), 
#       "hazard), Scale:", scales[class], "\n")
# }