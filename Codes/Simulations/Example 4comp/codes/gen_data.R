
library("lcmm")
library("NormPsy")
library(tictoc)
paquidS <- paquid[which(paquid$agedem> paquid$age_init),]
paquidS$normMMSE <-normMMSE(paquidS$MMSE)
paquidS$age65 <- (paquidS$age - 65)/10
seed <- 1; n <- 1000; censor_rate <- 0.05

# mj1 <- jlcmm(fixed = normMMSE ~ age65 +I(age65^2) + CEP,
#              random =~ age65+I(age65^2),
#              survival = Surv(age_init, agedem,dem)~CEP+male, hazard = "Weibull",
#              subject = 'ID', data = paquidS, verbose = F, logscale = T)
# summary(mj1)  


# tic()

# mj4b <- gridsearch(rep = 75, maxiter = 30, minit = mj1, cl = parallel::detectCores() - 1,
#                   jlcmm(fixed = normMMSE ~ age65 + I(age65^2) + CEP,
#                         random = ~ age65 + I(age65^2), 
#                         mixture = ~ age65 + I(age65^2), 
#                         ng = 4,
#                         survival = Surv(age_init, agedem, dem) ~ CEP + male, 
#                         hazard = "Weibull",
#                         subject = 'ID', data = paquidS,
#                         B = mj1, logscale = T)
#                   )

# toc()
# save(mj4b,file= "init4model.RData")
load("init4model.RData")
summary(mj4b)

init <- estimates(mj4b, cholesky=F)
D <- VarCovRE(mj4b)
#  matrix(D[,1], nrow=3, ncol=3, byrow = F)
cholesky <- metafor::vec2mat(D[,1], diag=T)
cholesky[3,1] <-cholesky[2,2] 
cholesky[2,2] <- cholesky[1,3]
cholesky[1,3] <- cholesky[3,1]
cholesky <- round(cholesky)
init[4:11] <- exp(init[4:11])
scales_est <- init[c(10,8,6,4)]^(-1/init[c(11,9,7,5)])
shapes_est <- init[c(11,9,7,5)]
betas_matrix <- round(matrix(init[14:25], ncol = 3))

residual_err <- 10

# Load required libraries
library(dplyr)
library(ggplot2)

# Set parameters for simulation
# set.seed(123)  # for reproducibility
# seed <- 1;n <- 500; censor_rate <- 0.05
set.seed(seed)
# source("NNEM/helper_functions/util.R")

# example 3: Xue & Yao 2022
n_samples <- n
x <- runif(n_samples, -5, 12)
# X <- cbind(1, x) 
x <- sort(x)
true_mixing_prop <- function(x) {
    pi1 <- (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4)) * as.numeric(x < 3)
    pi2 <- (1 - (2 * exp(-0.1 * x^4)) / (1 + exp(-0.1 * x^4))) * as.numeric(x < 0)
    pi3 <- (2 * exp(-0.1 * (x - 6)^4)) / (1 + exp(-0.1 * (x - 6)^4)) * as.numeric(x >= 3)
    pi4 <- 1 - (pi1 + pi2+ pi3)
    return(cbind(pi1, pi2, pi3, pi4))
  }
library(dplyr)
pi_init <- true_mixing_prop(x)

df <- data.frame(x,pi_init)
df<-df|>arrange(x)
plot(df$x,df$pi1, type="l", col="blue")
lines(df$x, df$pi2)
lines(df$x,df$pi3, col="green")
lines(df$x,df$pi4, col="red")
component_assignments <- apply(pi_init, 1, function(p) sample(1:4, 1, prob = p))
# component_assignments <- sample(1:3,size=n, prob=c(.3,.4, .3), replace = T)
# probs <- matrix(rep(c(.3,.4, .3), n), ncol = 3, byrow = T)
# pi_init <- probs
  shapes <- round(shapes_est)
  scales <- round(scales_est)
#  scales <- scales^(-1/shapes)
  # parms <- get_parms("weibulli")
  # Generate data
  times <- numeric(n)
  # i <- 1
  # lam_D <- parms$lam_D
  time_L <- runif(n, min=0, max=1)
  Nsub <- n
  m <- length(shapes)
  g <- components <- component_assignments; 
  table(components)/n
  # slopes <- parms$slopes
  library(flexsurv)
  # print(-1*slopes)
#   g <- components
  # ebx <- exp(rowSums(slopes[g,] * cbind(X3,X4, X5))) # proportional hazards term
  for (i in 1:m) {
    component_indices <- which(components == i)
    if (length(component_indices) > 0) {
      # Calculate linear predictor
      shape_g <- shapes[i]
      scale_g <- scales[i]
      # adjusted_scales <- scales[i] * exp(linear_pred[component_indices])
        #     times[component_indices] <- gen_model_survival(ebx[component_indices],
        # dist = "weibulli", 
        # scale=scale_g, 
        # shape=shape_g)
      times[component_indices] <- rweibull(length(component_indices), shape=shape_g, scale=scale_g)
    }
  }
  hist(times)
  # scales <- 
  # censor_rate <- "0.2"; #!!NOTE censor_rate=0.2 is more like 90% and the if you don't multiply by 10 and it's about 50% if you multply lambda by 10; censor_rate=0.5 yields an actual censor rate of 70% if you multiply by 10 and an actual censor rate of 0.25% if you don't multiply 
  # r <- 0.05 # 5%
  # # r <- 0.2 # 25%
  # # r <- 0.5 # 50%
  # tmp_parms <- parms
  time_L <- numeric(Nsub)

      if (censor_rate==0){ 
        time_C <- Inf
    } else if(censor_rate==0.05){
        time_C <- rweibull(n, shape =5, scale = 150)
    }else if(censor_rate==0.25){
      time_C <- rweibull(n, shape = 3, scale = 125)
    }else if(censor_rate==0.5){
      time_C <- rweibull(n, shape = 2, scale = 100)
    }
  delta <- as.numeric(time_C  >= times); print(table(delta)/n);
  time_Y <- pmin(times, time_C)
  hist(time_Y)
  # hist(time_C)

# Model parameters
n_subjects <- n     # number of individuals
# n_timepoints <-     # number of measure_time points per individual

hist(table(paquidS$ID))
median(table(paquidS$ID))
# Generate the data structure
data <- data.frame()
for(i in seq_len(n_subjects)){
  # i <- 1
n_timepoints <- rpois(1, 4)
data_i <- expand.grid(
  ID = i,
  # measure_time = seq(from=0,to=(n_timepoints), by=6)
  measure_time = sort(c(0.1,rweibull(n_timepoints, 2, scale=1.5)))

)
data <- rbind(data,data_i)}

# Generate individual random effects (random intercepts)
# random_effects <- data.frame(
#   ID = 1:n_subjects,
#   u_i = rnorm(n_subjects, mean = 0, sd = sigma_u)
# )

datalong <- data.frame()
for(i in 1:length(times)){

data_tmp <- data[data$ID==i,]
# data_tmp <- data_tmp[data_tmp$measure_time<times[i],]
data_tmp$event_time <- time_Y[i]
data_tmp$delta <- delta[i]
datalong <- rbind(datalong, data_tmp)
}
data <- datalong
sd_e <- residual_err
library(mvtnorm)
    mu_g <- matrix(0, 2, 2)
    rho <- 0.5
    cov <- rho*sqrt(9)*sqrt(4)
    D_g_og <- cholesky
    matrixcalc::is.positive.semi.definite(D_g_og) 
    # D_g_og <- diag(2)
    # prop_0 <- 1
    # prop_1 <- 0.5
    # prop_2 <- 1.5
    # prop_3 <- 2.5
    #? for the proportions, one must be very careful with the ordering particularly for initialization otherwise you will run into issues with estimation
    # props <- c(prop_3, prop_2, prop_0)
    # D_matrix2 <- (prop_2^2)*D_g_og
    # D_matrix3 <- (prop_3^2)*D_g_og
#  i <- 1


    
    ranef <- matrix(0, nrow=Nsub,ncol=3)
bi <- ranef
    datalong$measure_time2 <- datalong$measure_time^2
    X <- model.matrix(~1+ measure_time +(measure_time2), data=datalong)
    Z <- X
data_new <- time_Y

    for(i in 1:m){
    component_indices <- which(components == i)
    if (length(component_indices) > 0) {
        D <- D_g_og
        # ranef <- pi_init[,i]*
        ranef <- pi_init[,i]*rmvnorm(n, sigma = D)
        bi <- bi+ranef
    # init_ranef <- ranef + init_ranef
    }
    }
    

# 

    # # fixef <- c(0.5,2,1)
    # betas_matrix <- matrix(c(0.5,2, 
    #                         3, 2,
    #                         2, 1), ncol = 2, byrow=T)
    # # betas_matrix <- matrix(c(fixef,fixef), nrow = 3, ncol=2)
    pseudo_g <- g[datalong$ID]
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
D <- D_g_og
  V_blocks <- list()
    # omega <- props
mu <- c()
  for (i in 1:n) {
    # i <- 1
    idx <- which(data$ID == i)
    g <- component_assignments[i]
    Z_i <- Z[idx, , drop = FALSE]
    V_i <- Z_i %*% (D) %*% t(Z_i) + (sigma_e^2) * diag(length(idx))
    X_i <- X[idx, ,drop=FALSE]
    XBg <- X_i %*% betas_matrix[g, ]
    V_blocks[[i]] <- V_i
    mu <- c(mu, XBg)
  }
V <- Matrix::bdiag(V_blocks)
library(mvnfast)
y <- t(rmvn(1, mu, V))
datalong$y <- y
head(datalong,10)
# Clean up intermediate columns (optional): don't need to do this if we generate from a poison
# data <- data %>%
#   select(ID, measure_time, y) %>%
#   arrange(ID, measure_time)

# Display first few rows
head(datalong, 15)
datalong$g <-pseudo_g 
# Summary statistics
cat("Summary of generated data:\n")
cat("Number of subjects:", n_subjects, "\n")
cat("Number of observations per subject:", mean(as.numeric(n_timepoints)), "\n")
cat("Total observations:", nrow(data), "\n")
cat("Mean response:", round(mean(datalong$y), 2), "\n")
cat("SD of response:", round(sd(datalong$y), 2), "\n")

# Plot a sample of trajectories
sample_subjects <- sample(1:n_subjects, min(50, n_subjects))
data_sample <- datalong %>% filter(ID %in% sample_subjects)

p <- ggplot(data_sample, aes(x = measure_time, y = y, group = ID, color = factor(ID))) +
  geom_line() +
  geom_point() +
  labs(title = "Sample Individual Trajectories",
       x = "Time",
       y = "Response (Y)",
       color = "Subject ID") +
  theme_minimal()
print(p)
# Verify the model by fitting it back to the generated data
library(lme4)
model_fit <- lmer(y ~ measure_time + (measure_time+measure_time2| ID), data = datalong)
print(summary(model_fit))


  data <- datalong

    m <- 4
    D_g <- D_g_og
    # alphas <- props
    # alphas
    # gammas_matrix <- slopes[1:m,]
    # mu_g <- mu_g
    nis<-table(data$ID)


library('data.table')
sim <- data.table(data)
temp2 <- subset(sim, select = c(ID, event_time, delta)) # baseline
data_wide <- unique(temp2)

X <- model.matrix(~1+ measure_time +(measure_time2), data=datalong)
Z <- X
bi <- list()
# nis <- table(num_measure)
#? compute class specific random effects 
m <- length(shapes)
# alphas
for(j in 1:m){
#   j <- 1
    sigma_g <- sd_e # error variances 

    D_g <- D_g_og
     # random effect variance-covariance matrices
    beta_g <- betas_matrix[j,] # fixed effects

  # i<-1
  bigs<-matrix(nrow = n, ncol = 3)
    for(i in 1:n){
      # i<-2
      n_i <- nis[i]
      # k <- 1
      sim <- data
      subset<-sim[sim$ID==i,]
      Z_i <- Z[sim$ID==i, , drop = FALSE]
      X_i <- Z_i
      p<-ncol(X_i)
      Xi_beta <- X_i %*% beta_g
      # xbs[i]<-list(Xi_beta)
      q<-ncol(Z_i)
      Y_i <- subset$y
      # ys<-list(Y_i)
      sigma<-sigma_g

      # 3by2 %*% 
      Vi <- as.matrix(Z_i%*% D_g%*% t(Z_i) + sigma_g^2 * diag(n_i))
      big <- D_g %*% t(Z_i) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
      bigs[i,] <- big
      }
  bi[j] <- list(bigs)
}



library(nnet)
cholesky <- t(chol(D_g_og))
tcrossprod(cholesky)
# omegas <- alphas[1:2]
theta_orig <- theta <- list(pj=pi_init, scale=scales, shape=shapes, 
   betas=betas_matrix, ranef=bi, 
  D=D_g_og, cholesky=cholesky, sigma_e=sd_e)
save(theta_orig, file="theta_orig_k4.RData")


# Plot true class-specific survival probabilities
library(ggplot2)
library(dplyr)

# Use your existing parameters
shapes <- round(shapes_est)
scales <- round(scales_est)

print("True Weibull parameters:")
print(data.frame(Class = 1:4, Shape = shapes, Scale = scales))

# Create time sequence for plotting
time_seq <- seq(0, 200, by = 1)

# Function to calculate Weibull survival probability
weibull_survival <- function(t, shape, scale) {
  exp(-(t/scale)^shape)
}

# Calculate survival probabilities for each class
survival_data <- data.frame()

for(class in 1:4) {
  class_survival <- data.frame(
    time = time_seq,
    survival_prob = weibull_survival(time_seq, shapes[class], scales[class]),
    class = paste("Class", class)
  )
  survival_data <- rbind(survival_data, class_survival)
}

# Create the plot
p1 <- ggplot(survival_data, aes(x = time, y = survival_prob, color = class)) +
  geom_line(size = 1.2) +
  labs(
    title = "True Class-Specific Survival Probabilities",
    subtitle = paste("Weibull distributions with shapes:", paste(shapes, collapse = ", "), 
                     "and scales:", paste(scales, collapse = ", ")),
    x = "Time",
    y = "Survival Probability",
    color = "Latent Class"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(65, 120), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c("Class 1" = "#1f77b4", "Class 2" = "#ff7f0e", 
                               "Class 3" = "#2ca02c", "Class 4" = "#d62728")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p1)

# Alternative: Base R plot
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
plot(time_seq, weibull_survival(time_seq, shapes[1], scales[1]), 
     type = "l", col = "blue", lwd = 2,
     ylim = c(0, 1), xlim=c(65, 120),xlab = "Time", ylab = "Survival Probability",
     main = "True Class-Specific Survival Probabilities")

for(class in 2:4) {
  lines(time_seq, weibull_survival(time_seq, shapes[class], scales[class]), 
        col = c("orange", "green", "red")[class-1], lwd = 2)
}

legend("topright", 
       legend = paste("Class", 1:4, "(λ=", scales, ", κ=", shapes, ")"),
       col = c("blue", "orange", "green", "red"),
       lwd = 2, cex = 0.8)

# Additional analysis: median survival times
median_survival <- numeric(4)
for(class in 1:4) {
  # Median survival time for Weibull: scale * (ln(2))^(1/shape)
  median_survival[class] <- scales[class] * (log(2))^(1/shapes[class])
}

print("Median survival times by class:")
print(data.frame(Class = 1:4, 
                 Shape = shapes, 
                 Scale = scales, 
                 Median_Survival = round(median_survival, 2)))

# Hazard ratios (relative to Class 1)
cat("\nHazard characteristics:\n")
for(class in 1:4) {
  cat("Class", class, "- Shape:", shapes[class], 
      "(", ifelse(shapes[class] > 1, "increasing", 
                  ifelse(shapes[class] == 1, "constant", "decreasing")), 
      "hazard), Scale:", scales[class], "\n")
}