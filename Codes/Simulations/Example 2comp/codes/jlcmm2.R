# Load required libraries
library(dplyr)
library(ggplot2)

# Set parameters for simulation
# set.seed(123)  # for reproducibility
seed <- 1;n <- 100; censor_rate<-0.05
set.seed(seed)
source("NNEM/helper_functions/util.R")

# example 3: Xue & Yao 2022
n_samples <- n
  u <- runif(n_samples, -7.5, 7.5)
  
  # Generate x_{ij} for j=1,...,9
  # x_1_9 <- matrix(runif(n_samples * 9, -1.5, 1.5), ncol = 9)
  
  # Generate x_{i10}
  # x_10 <- u - 0.5 * rowSums(x_1_9)
  
  # Full design matrix
  # X <- cbind(1, x_1_9, x_10)  # (n_samples, 11)
  
  # True beta coefficients (length 11)
  # beta1 <- c(0, 2, 1, 3, 1, 2, 3, 2, 1, 1, 1.5)
  # beta2 <- c(0, 1, 3, 3, 4, 6, 3, 4, 6, 5, 3)
    x=sort(runif(n,0,1))
    k=2
    mix.prop=0.1+0.8*sin(pi*x);mix.prop=cbind(1-mix.prop,mix.prop) #Huang&Yao(2012)#2
    #mix.prop=0.7+0.3*sin(pi*x);mix.prop=cbind(1-mix.prop,mix.prop) #Huang&Yao(2012)#3
    #mix.prop=1-((x-0.5)/1)^2;mix.prop=cbind(1-mix.prop,mix.prop) #Young&Hunter(2010)
    #mix.prop=2*exp(-0.1*x^4)/(1+exp(-0.1*x^4));mix.prop=cbind(1-mix.prop,mix.prop) #Xue&Yao(2022)
    #mix.prop=1/(1+exp(-0.5*x));mix.prop=cbind(1-mix.prop,mix.prop) #Young&Hunter(2010)
    mix.sigma2 = c(0.09,0.16)
    mix.mu = cbind(4-2*x,3*x)
    mix.beta = cbind(c(4,-2),c(0,3))
    z = sapply(1:n,function(i) sample(1:k,1,replace=T,prob=mix.prop[i,]))
    component_assignments <- z
    pi_init <- mix.prop
df <- data.frame(x,pi_init)
# df<-df|>arrange(x)
plot(df$x,df$mix.prop, type="l", col="blue")
lines(df$x,df$V1, col="green")
lines(df$x,df$pi4, col="red")
# component_assignments <- apply(pi_init, 1, function(p) sample(1:3, 1, prob = p))
# component_assignments <- sample(1:3,size=n, prob=c(.3,.4, .3), replace = T)
# probs <- matrix(rep(c(.3,.4, .3), n), ncol = 3, byrow = T)
# pi_init <- probs
  shapes <- c(2, 1.01)
  scales <- c(4.5, 50)
#  scales <- scales^(-1/shapes)
  parms <- get_parms("weibulli")
  # Generate data
  times <- numeric(n)
  # i <- 1
  lam_D <- parms$lam_D
  time_L <- runif(n, min=0, max=1)
  Nsub <- n
  m <- length(shapes)
  g <- components <- component_assignments; 
  table(components)/n
  slopes <- parms$slopes
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
#  censor_rate <- 0.5 #!!NOTE censor_rate=0.2 is more like 90% and the if you don't multiply by 10 and it's about 50% if you multply lambda by 10; censor_rate=0.5 yields an actual censor rate of 70% if you multiply by 10 and an actual censor rate of 0.25% if you don't multiply 
#  r <- 0.05; #5%
# r <- 0.55; #25%
# # censor_rate <- "0.5"; r<-0.6; #50%
# r <- 1 # 50%
target_censor_rate <- 0.1
  tmp_parms <- parms
  time_L <- numeric(Nsub)

      if (censor_rate==0){ 
        time_C <- Inf
    } else if(censor_rate==0.05){
        time_C <- rweibull(n, shape =10, scale = 150)
    }else if(censor_rate==0.25){
      time_C <- rweibull(n, shape = 10, scale = 50)
    }else if(censor_rate==0.5){
      time_C <- rweibull(n, shape = 0.5, scale = 25)
    }

  delta <- as.numeric(time_C  >= times); print(table(delta)/n);
  time_Y <- pmin(times, time_C)
  hist(time_Y)
  # hist((time_C))

# Model parameters
n_subjects <- n     # number of individuals
# n_timepoints <-     # number of measure_time points per individual
beta_0 <- 10          # fixed intercept
beta_1 <- 2           # fixed slope for measure_time
sigma_u <- 2          # standard deviation of random intercepts
# sigma_e <- 1          # standard deviation of residual errors

# Generate the data structure
data <- data.frame()
tot <- c(0,1,3,6,12,18,24) # include a baseline recording at 0
for(i in seq_len(n_subjects)){
  # i <- 1
n_timepoints <- rpois(1, 8)
  if(n_timepoints<=length(tot)){
  measure_times <- tot[1:n_timepoints]}
  else{
    measure_times <- tot
  }
data_i <- expand.grid(
  ID = i,
  # measure_time = seq(from=0,to=(n_timepoints), by=6)
  # measure_time =runif()

  measure_time = measure_times/12

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
data_tmp <- data_tmp[data_tmp$measure_time<times[i],]
data_tmp$event_time <- time_Y[i]
data_tmp$delta <- delta[i]
datalong <- rbind(datalong, data_tmp)
}

    sd_e <- sqrt(60)
    library(mvtnorm)
    mu_g <- matrix(0, 2, 2)
    rho <- 0
    cov <- rho*sqrt(9)*sqrt(4)
    # D_g_og <- matrix(c(9, cov, cov, 4), ncol = 2, nrow=2) 
    D_g1 <- 50
    # D_g2 <- 
data_new <- times

    # for(i in 1:m){
    # component_indices <- which(components == i)
    # if (length(component_indices) > 0) {
    #     D <- (props[i]^2)*D_g_og
    #     # ranef <- pi_init[,i]*
    #     ranef[component_indices,] <- rmvnorm(length(component_indices), sigma = D)
    #     # bi <- bi+ranef
    # # init_ranef <- ranef + init_ranef
    # }
    # }
    

# 

    # # fixef <- c(0.5,2,1)
    betas_matrix <- matrix(c(170,88, 
                            100, 1.2), ncol = 2, byrow=T)
    # # betas_matrix <- matrix(c(fixef,fixef), nrow = 3, ncol=2)
    pseudo_g <- g[datalong$ID]
    # # ranefsx <- bi[[pseudo_g]]
    # ranefs <- ranef[datalong$ID,]
    
    D <- matrix(D_g1)
    # Zb <- ranefs[,1] + ranefs[,2]*datalong$measure_time
    # y <- betas_matrix[,1][pseudo_g]+betas_matrix[,2][pseudo_g]*datalong$measure_time+ Zb + rnorm(nrow(datalong), sd=sd_e) 
library(MASS)
library(matrixcalc)
n_timepoints <- table(datalong$ID)
# Create design matrix
Z <- model.matrix(~1, data = datalong)
X <- model.matrix(~measure_time,data=datalong)
V_blocks <- list()
  #  omega <- props
mu <- c()
  for (i in 1:n) {
    # i <- 1
    idx <- which(datalong$ID == i)
    g <- component_assignments[i]
    Z_i <- Z[idx, , drop = FALSE]
    V_i <- ((Z_i %*% (D) %*% t(Z_i)) + ((sd_e^2) * diag(length(idx))))
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
cat("Average number of observations per subject:", mean(n_timepoints), "\n")
cat("Total observations:", nrow(data), "\n")
cat("Mean response:", round(mean(datalong$y), 2), "\n")
cat("SD of response:", round(sd(datalong$y), 2), "\n")

# Plot a sample of trajectories
sample_subjects <- sample(1:n_subjects, min(30, n_subjects))
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
model_fit <- lmer(y ~ measure_time + (1| ID), data = datalong)
print(summary(model_fit))


  data <- datalong

    m <- length(D)
    D_g <- D
    # alphas <- props
    # alphas
    # gammas_matrix <- slopes[1:m,]
    mu_g <- mu_g
    nis<-table(datalong$ID)


library('data.table')
sim <- data.table(data)
temp2 <- subset(sim, select = c(ID, event_time, delta)) # baseline
data_wide <- unique(temp2)
q <- ncol(Z)
bi <- list()
# nis <- table(num_measure)
#? compute class specific random effects 
m <- length(shapes)
# alphas
for(j in 1:m){
  #  j <- 1
    sigma_g <- sd_e # error variances 

    D_g <- D
     # random effect variance-covariance matrices
    beta_g <- betas_matrix[j,] # fixed effects

  #  i<-1
  bigs<-matrix(nrow = n, ncol = q)
    for(i in 1:n){
    #   i<-2
      n_i <- nis[i]
      # k <- 1
      sim <- data
      subset<-sim[sim$ID==i,]
      Z_i <- Z[sim$ID==i, , drop = FALSE]
      X_i <- X[sim$ID==i, ,drop=FALSE]
       p<-ncol(X_i)
      Xi_beta <- X_i %*% beta_g
      # xbs[i]<-list(Xi_beta)
      # Z_i<-
      q<-ncol(Z_i)
      Y_i <- subset$y
      # ys<-list(Y_i)
      sigma<-sigma_g

      # 3by2 %*% 
      Vi <- as.matrix(Z_i%*% D_g%*% t(Z_i) + sigma_g^2 * diag(n_i))
      big <- D_g %*% t(Z_i) %*% matrixcalc::matrix.inverse(Vi) %*% (Y_i - Xi_beta)
      bigs[i] <- big
      }
  bi[j] <- list(bigs)
}
# m
pi_init <- matrix(1/m, nrow = length(x), ncol = m)

print(pi_init)
library(nnet)
cholesky <- t(chol(D))
# tcrossprod(cholesky)
# omegas <- alphas[1:2]
theta <- list(pj=pi_init, scale=scales, shape=shapes, 
  betas=betas_matrix, ranef=bi, cholesky=cholesky,
  D=D, sigma_e=sd_e)
