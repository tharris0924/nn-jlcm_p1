  source("NNEM/Application/gen_data.R")
  source("NNEM/Application/em_utils.R")
  source("NNEM/Application/estep.R")
  source("NNEM/Application/mstep.R")

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
# gammas_params <- gammas_matrix
get_pred_surv <- function(shape_params, scale_params, test_data, mix.prob, survform, gammas_params) {
    survform <- update(survform, ~.-1)

 #* determine class belongings
  test_data <- data.table(test_data)
  test_data$rownames <- c(1:nrow(test_data))
   rowtokeep <- test_data[, rownames[1], by = ID][, V1]
   test_data <- test_data[rowtokeep,]
  # colnames(mix.prob) <- c("V1","V2")
    # pred <- matrix(0, nrow = N, ncol = k)
  library(tidyverse)
  # mix.prob <- pij[[1]]
  dat <- cbind(mix.prob, test_data)
   dat <- dat %>% arrange((event_time))
    k <- length(as.matrix(shape_params))
     ti <- c(dat$event_time)
  brier_surv <- Surv(time=ti, event=dat$delta)
    W <- model.matrix(survform, dat)

  dat <- as.matrix(dat)
   mix.prob <- as.matrix(dat[,c(1:k)])
  #  ti <- seq(min(test_data$event_time), max(test_data$event_time), length.out = 101)
  # mbvars
  # classmb <- paste0("~",mbvars, sep="")
  
      N <- nrow(test_data)

  # initial density function for Weibull
  pred <- matrix(0,nrow=length(ti), ncol=k )


      # ti <- test_data$event_time
      # status <- result$data_wide$delta

  Shat <- NULL

    for(i in 1:N){    
        pred_tmp <- matrix(0, nrow=k, ncol=length(ti))
      for (j in 1:k) {
    # j <- 1
    shape <- shape_params[j]
    scale <- scale_params[j]
    gammas <- gammas_params
    linear_pred <- exp((W[i,]%*%gammas)/shape) #updated from Prof. Breheny's lecture slides p18 (Weibull distribution)
    # linear_pred <- exp()
    # delta <- delta
    new_scales <- (linear_pred)*scale
    # fixed effects in survival component
    # wxi <- as.matrix(w)
    pred_tmp[j,]<-pweibull(ti,shape, new_scales, lower.tail = F)
    }
  Shat_raw <- colSums(pred_tmp * mix.prob[i,])
  #  Shat_raw <- sort(Shat_raw, decreasing = T)
    Shat <- rbind(Shat, Shat_raw)
    }
  return(list(Shat=Shat, ti=ti, brier_surv=brier_surv))
}

# plot(ti,IBS$bsc, type="l")

library(pec)

pred_yis <- function(data_long, theta, long_form_g, long_form_o,ranform) {

  betas_matrix <- theta$betas
  sigmas <- theta$sigma_e
  # D_matrices <- theta$D
  beta_o <- theta$beta_o # overall sample fixed effects parameters

  # data prep
  y<- model.frame(long_form_g,data_long)[,1]
  Xg <- model.matrix(long_form_g, data_long)
  Xo <- model.matrix(long_form_o, data_long)
  Zg <- model.matrix(ranform, data_long)
  data_long$y <- y

  N <- length((data_long$ID)) #;print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix) #;print(m);print(betas_matrix) # number of unique groups for mixture model
  pred_yi <- matrix(0,nrow = N, ncol = m);
  nis<-table(data_long$ID)
  IDs<-unique(data_long$ID)
  # length(nis)
  D_matrices <- tcrossprod(theta$cholesky);# print(D_matrices)
  for (k in 1:m) {
    # k <- 1
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects
   pred_yij <- NULL
      bi <- theta$ranef[[k]]

    
    for(i in seq_len(length(IDs))){

      # i<-2
      # n_i <- nis[i]
      # k <- 1
      ids <- IDs[i] 
      # data_long <- data;i <- 1
      subset <- data_long[data_long$ID == ids, ] #! ensure that IDs are consistently named across datasets
            # subset<-sim[sim$ID==i,]
      Z_ig <- Zg[data_long$ID == ids, , drop = FALSE]    # Random Effects Design Matrix
      X_ig <- Xg[data_long$ID==ids,,drop=FALSE ]# Cluster-specific Fixed Effects Design Matrix
      X_io <- Xo[data_long$ID==ids,,drop=FALSE] # Overall Fixed Effects Design Matrix (*no intercept)
      p <- ncol(X_ig)
      Xi_beta_k <- matrix(X_ig %*% beta_g + X_io %*% beta_o + Z_ig %*%bi[i,])
      pred_yij <- rbind(pred_yij,Xi_beta_k)
    }
    pred_yi[,k] <- pred_yij
  }
  return(pred_yi)
}
# ranef<-pred_ranef(data_long=test_data_long,theta=result$theta,long_form_g,long_form_o,ranform)
pred_ranef<-function(data_long, theta, long_form_g, long_form_o,ranform) {

  betas_matrix <- theta$betas
  sigmas <- theta$sigma_e
  # D_matrices <- theta$D
  beta_o <- theta$beta_o # overall sample fixed effects parameters

  # data prep
  y<- model.frame(long_form_g,data_long)[,1]
  Xg <- model.matrix(long_form_g, data_long)
  Xo <- model.matrix(long_form_o, data_long)
  Zg <- model.matrix(ranform, data_long)
  data_long$y <- y
  q <- ncol(Zg)
  N <- length((data_long$ID)) #;print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix) #;print(m);print(betas_matrix) # number of unique groups for mixture model
  pred_yi <- matrix(0,nrow = N, ncol = m);
  n <- length(unique(data_long$ID))
  nis<-table(data_long$ID)
  IDs<-(unique(data_long$ID))
  bi <- list()
  # length(nis)
  D_matrices <- tcrossprod(theta$cholesky);# print(D_matrices)
  for (k in 1:m) {
    #  k <- 1
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects
      # bi <- theta$ranef[[k]]

        bigs<-matrix(nrow = n, ncol = q)

    for(i in seq_len(length(IDs))){

      #  i<-2
      # n_i <- nis[i]
      # k <- 1
      ids <- IDs[i] 
      # data_long <- data;i <- 1
      subset <- data_long[data_long$ID == ids, ] #! ensure that IDs are consistently named across datasets
            # subset<-sim[sim$ID==i,]
      Z_ig <- Zg[data_long$ID == ids, , drop = FALSE]    # Random Effects Design Matrix
      X_ig <- Xg[data_long$ID==ids,,drop=FALSE ]# Cluster-specific Fixed Effects Design Matrix
      X_io <- Xo[data_long$ID==ids,,drop=FALSE] # Overall Fixed Effects Design Matrix (*no intercept)
      p <- ncol(X_ig)
      Xi_beta <- matrix(X_ig %*% beta_g + X_io %*% beta_o)
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
  bi[k] <- list(bigs)
  }
 
  return(bi)
}



# library(JM)
# ?aucJM

# # install.packages('nftbart')
# library(pec)



# plot(survfit(itx~1))

plot_survival_curves <- function(model = NULL, data_wide, mix_prob = NULL, 
                                  shape_params = NULL, scale_params = NULL, gammas_params = NULL,
                                  survform = NULL, title = "Survival Curves", 
                                  xlim = c(65, 100), ylim = c(0, 1)) {
  
  
  # Load required libraries
  library(tidyverse)
  library(survival)
  library(broom)
  library(ggplot2)
  
  # Extract parameters based on input type
  if (!is.null(model)) {
    # Extract parameters from JLCMM model
    sample_init <- estimates(model, cholesky = FALSE)
    
    # Extract Weibull parameters
    idWeib_shape <- which(grepl('Weibull2', names(sample_init)))
    idWeib_scale <- which(grepl('Weibull1', names(sample_init)))
    scale_params <- exp(sample_init[idWeib_scale])^(-1/exp(sample_init[idWeib_shape]))
    shape_params <- exp(sample_init[idWeib_shape])
    gammas_params <- sample_init[(max(idWeib_shape)+1):(max(idWeib_shape)+2)]
    
    # Get mixing probabilities from model if not provided
    if (is.null(mix_prob)) {
      mix_prob <- model$pprob[, -c(1:2)]  # Remove ID and class columns
    }
  } else {
    # Use provided parameters
    if (is.null(shape_params) || is.null(scale_params) || is.null(gammas_params)) {
      stop("Either provide a model or all parameter vectors (shape_params, scale_params, gammas_params)")
    }
    if (is.null(mix_prob)) {
      stop("mix_prob must be provided when using custom parameters")
    }
  }
  survform <- update(survform, ~.-1)
  N <- nrow(data_wide)
  k <- nrow(as.matrix(shape_params))

  # f_td <- matrix(0,nrow = N, ncol = k)
  # w <- model.matrix(~ -1+X3+X4+X5, data=data_wide)
  delta <- data_wide$delta
  # initial density function for Weibull
  pred <- matrix(0, nrow = N, ncol = k)
  dat <- cbind(mix_prob,data_wide)
  # colnames(dat) <- c("ID", "event_time", "delta", "V1", "V2")
   dat <- dat %>% arrange((event_time))
   mix.prob <- as.matrix(dat[,1:k])
  W <- model.matrix(survform, dat)

  # ti <- seq(0, max(data_wide$event_time), length.out = 100)
  ti <- dat$event_time


  pred <- matrix(0, nrow = N, ncol = k)
    for (j in 1:k) {
    # j <- 1
    shape <- shape_params[j]
    scale <- scale_params[j]
    gammas <- gammas_params
    linear_pred <- exp((W%*%gammas)/shape) #updated from Prof. Breheny's lecture slides p18
    # linear_pred <- exp()
    delta <- delta
    new_scales <- (linear_pred)*scale
    # fixed effects in survival component
    # wxi <- as.matrix(w)
    for(i in 1:N){
    pred[i,j]<-pweibull(ti[i],shape, new_scales[i], lower.tail = F)
  }}
  
  pred[,1:k] <- sapply(1:k, function(i) sort(pred[,i], decreasing = TRUE))  

  
  # Prepare predicted curves data for ggplot
  pred_data <- data.frame(
    ti = rep(ti, k),
    survival = c(pred),
    cluster = rep(paste("Cluster", 1:k), each = length(ti)),
    type = "Predicted"
  )
  
  # Prepare Kaplan-Meier data with confidence intervals
  km_data <- data.frame()
  
  for (i in 1:k) {
    km_fit <- survfit(Surv(dat$event_time, dat$delta) ~ 1, weights = mix.prob[, i])
    km_tidy <- tidy(km_fit)
    km_tidy$cluster <- paste("Cluster", i)
    km_tidy$type <- "Kaplan-Meier"
    km_data <- rbind(km_data, km_tidy)
  }
  
  # Create the plot
  p <- ggplot() +
    # Predicted curves
    geom_line(data = pred_data, aes(x = ti, y = survival, color = cluster), 
              linetype = "solid", size = 1) +
    # Kaplan-Meier confidence intervals
    geom_ribbon(data = km_data, 
                aes(x = time, ymin = conf.low, ymax = conf.high, fill = cluster), 
                alpha = 0.2) +
    # Kaplan-Meier curves  
    geom_line(data = km_data, 
              aes(x = time, y = estimate, color = cluster), 
              linetype = "dashed", size = 1) +
    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +
    labs(x = "Age in Years", y = "Survival Probability", 
         color = "Cluster", fill = "Cluster") +
    theme_minimal() +
    ggtitle(title)
  
  # Return both the plot and the underlying data
  return(list(
    plot = p,
    predicted_data = pred_data,
    kaplan_meier_data = km_data,
    group_spec_surv=pred,
    parameters = list(
      shape_params = shape_params,
      scale_params = scale_params,
      gammas_params = gammas_params
    )
  ))
}


library(patchwork)
# combined_plot_patchwork <- p1/ p2 +
#   plot_annotation(
#     title = "Class-Specific Survival Probabilities",
#     theme = theme(plot.title = element_text(size = 16, face = "bold"))
#   )
# combined_plot_patchwork

library(SurvMetrics)
#! NNEM EVAL
load("NNEM/Application/nnjlcmk2.RData")
ll1 <- result$loglik
k <- 20
(BIC2=-2*ll1+k*log(nrow(result$data_wide)))
(AIC2=-2*ll1+2*k)
result$theta_sd$`Longitudinal Submodel`

shape_params <- result$theta$shape
scale_params <- result$theta$scale
gammas_params <- result$theta$gamma_matrix
result2 <- plot_survival_curves(
  data_wide = result$data_wide,
  mix_prob = result$postprob,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "NN-JLCM"
)
print(result2$plot)
mix.prob <- result$postprob
stat_i <- matrix(ncol=3, nrow=nrow(result2$group_spec_surv))
cumHaz <- -log(result2$group_spec_surv)
stat <- 0
  for(i in 1:nrow(result$data_wide)){
    sums <- 0
    for(j in 1:ncol(result$theta$pj)){
    # i <- j <- 1
    big <- result$theta$ranef[[j]][i,]
    sum <- mix.prob[i,j]*(result$data_wide$delta[i]-cumHaz[i,j])*big
    sums <- sum+sums
    # stat <- sum+stat
    }
  stat_i[i,] <- sums
  }
# mean_stat <- tcrossprod(stat)/n
n <- nrow(result$data_wide)
var <- matrix(0, nrow=3, ncol=3)
for(i in 1:nrow(result$data_wide)){
    var <- var + (tcrossprod(stat_i[i,]))
}
var <- var-tcrossprod(colSums(stat_i))/n
(score_stat2 <- t(colSums(stat_i))%*%matrixcalc::matrix.inverse(var)%*%colSums(stat_i))
(scoretest_pvalue2 <- 1-pchisq(score_stat2, 3))
#! training predictions
# estep.result <- estep(prior_mixing=result$theta$pj, result$theta, long_form_o, long_form_g, ranform, survform, result$data_long, result$data_wide)
mix.prob <- result$postprob
entropy <- sum(mix.prob*log(mix.prob))
(ICL1_k2 <- BIC2 - entropy)
(ICL2_k2 <- BIC2 -2*sum(log((apply(mix.prob, MARGIN = 1, FUN = max)))))
measure_entropy_2 <- 1 + (sum(entropy))/(nrow(result$data_wide)*log(ncol(mix.prob)))
    classif <- NULL
    classifY <- NULL
    cl.table <- NULL
    thr.table <- NULL
threshold=c(0.7,0.8,0.9)
pprob <- cbind(result$C,mix.prob)
ng <- ncol(mix.prob)
    for (g in 1:ng) 
    {
    temp<- subset(pprob,pprob[,1]==g)
    # tempY<- subset(pprobY,pprobY[,1]==g)
    temp1<-apply(temp[,1+1:ng],2,mean)
    # temp1Y<-apply(tempY[,1+1:ng],2,mean)
    cl.table<-rbind(cl.table,temp1)
    classif<-cbind(classif,length(temp[,1]))
    # classifY<-cbind(classifY,length(tempY[,1]))
    if(!is.null(threshold)) thr.table <- cbind(thr.table,sapply(threshold,function(x) length(which(temp[,1+g]>x)))/length(temp[,1]))  
   }
   
   classif <- rbind(classif,100*classif/nrow(result$data_wide))
  #  classifY <- rbind(classifY,100*classifY/x$ns)
      if(!is.null(thr.table))
   {
    thr.table <- 100*thr.table  
    
    rownames(thr.table) <- paste("prob>",threshold,sep="") 
    colnames(thr.table) <- paste("class",1:ng,sep="") 
   }
    if(!is.null(thr.table))
    {
     cat("Posterior probabilities above a threshold (%):", "\n")
     print(round(thr.table,2))
     cat(" \n")    
    }
   rownames(cl.table)<-paste("class",1:ng,sep="")
   colnames(cl.table)<-paste("prob",1:ng,sep="")
   colnames(classif)<-paste("class",1:ng,sep="")
   rownames(classif)<-c("N","%")
  #  colnames(classifY)<-paste("class",1:ng,sep="")
  #  rownames(classifY)<-c("N","%")
    cat(" \n")
    cat("Posterior classification based on longitudinal and time-to-event data:", "\n")
    print(round(classif,2))
    cat(" \n")
    
    cat("Posterior classification table:", "\n")
    cat("     --> mean of posterior probabilities in each class", "\n")
    print(round(cl.table,4))

cl.2 <- cl.table
# mix.prob <- estep.result$P
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = result$data_wide, mix.prob = mix.prob, survform = survform)
library(SurvMetrics)
it2 <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)

rownames(mix.prob) <- result$data_wide$ID
yhat <- pred_yis(data_long=result$data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
preds <- rowSums(yhat * mix.prob[as.character(result$data_long$ID),])
local_MSE_y_in <- Metrics::mse(actual = result$data_long$y, predicted = preds)


#! testing predictions
  n <- nrow(test_data_wide)
  # m <- nrow(betas_matrix)
  x <- model.frame(gform, result$data_wide)
# mixing proportions using nnet

    # For more than 2 components, use multinomial
    # postprob <- e_step_result$postprob

    nn_model <- nnet::nnet(
      x = x,          # Input features
      y = mix.prob,          # Target is component probabilities
      size = 10,      # Hidden layer size
      decay = 0.01,   #! Weight decay: this is something which can be tweaked
      maxit = 5000,   # Max iterations
      trace = FALSE,  # No output during training
      softmax = TRUE,  # Use softmax for multinomial
      entropy=TRUE # maximum conditional likelihood fitting
    )
    
    # Predict probabilities (this code doesn't work for some stupid reason, 
  # extract directly from the nnet object)
x_new <-  model.frame(gform, test_data_wide)
    prob <- predict(nn_model, x_new)
    # new_prior <- fitted(nn_model)
estep.result <- estep(prior_mixing=prob, theta=result$theta, long_form_o, long_form_g, ranform, survform,data_long =  test_data_long,data_wide =  test_data_wide)
mix.prob <- estep.result$P
# mix.prob <- estep.result$P
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = test_data_wide, mix.prob = mix.prob, survform = survform)

it2_test <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)


rownames(mix.prob) <- test_data_wide$ID

ranef <- pred_ranef(data_long=test_data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
result$theta$ranef <- ranef
yhat <- pred_yis(data_long=test_data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
preds <- rowSums(yhat * mix.prob[as.character(test_data_long$ID),])
local_MSE_y_out_2 <- Metrics::mse(actual = test_data_long$y, predicted = preds)


load("NNEM/Application/nnjlcmk3.RData")
k <- nrow(result$theta_sd$`Longitudinal Submodel`)+nrow(result$theta_sd$`Survival Submodel`)
rbind(result$theta_sd$`Longitudinal Submodel`, result$theta_sd$`Survival Submodel`)
ll3 <- result$loglik
(BIC3=-2*ll3+k*log(nrow(result$data_wide)))
(AIC3=-2*ll3+2*k)

shape_params <- result$theta$shape
scale_params <- result$theta$scale
gammas_params <- result$theta$gamma_matrix
#! training predictions
# estep.result <- estep(prior_mixing=result$theta$pj, result$theta, long_form_o, long_form_g, ranform, survform, result$data_long, result$data_wide)
mix.prob <- result$postprob
# mix.prob <- estep.result$P

# Assuming you have extracted parameters from your other model:
result2 <- plot_survival_curves(
  data_wide = data_wide,
  mix_prob = result$postprob,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "In-Sample NN-JLCM"
)
nnjm_train <- result2

print(result2$plot)
pprob <- cbind(result$C, mix.prob)
    classif <- NULL
    classifY <- NULL
    cl.table <- NULL
    thr.table <- NULL
threshold=c(0.7,0.8,0.9)
ng <- ncol(mix.prob)
    for (g in 1:ng) 
    {
    temp<- subset(pprob,pprob[,1]==g)
    # tempY<- subset(pprobY,pprobY[,1]==g)
    temp1<-apply(temp[,1+1:ng],2,mean)
    # temp1Y<-apply(tempY[,1+1:ng],2,mean)
    cl.table<-rbind(cl.table,temp1)
    classif<-cbind(classif,length(temp[,1]))
    # classifY<-cbind(classifY,length(tempY[,1]))
    if(!is.null(threshold)) thr.table <- cbind(thr.table,sapply(threshold,function(x) length(which(temp[,1+g]>x)))/length(temp[,1]))  
   }
   
   classif <- rbind(classif,100*classif/nrow(result$data_wide))
  #  classifY <- rbind(classifY,100*classifY/x$ns)
      if(!is.null(thr.table))
   {
    thr.table <- 100*thr.table  
    
    rownames(thr.table) <- paste("prob>",threshold,sep="") 
    colnames(thr.table) <- paste("class",1:ng,sep="") 
   }
    if(!is.null(thr.table))
    {
     cat("Posterior probabilities above a threshold (%):", "\n")
     print(round(thr.table,2))
     cat(" \n")    
    }
    
   rownames(cl.table)<-paste("class",1:ng,sep="")
   colnames(cl.table)<-paste("prob",1:ng,sep="")
   colnames(classif)<-paste("class",1:ng,sep="")
   rownames(classif)<-c("N","%")
  #  colnames(classifY)<-paste("class",1:ng,sep="")
  #  rownames(classifY)<-c("N","%")
    cat(" \n")
    cat("Posterior classification based on longitudinal and time-to-event data:", "\n")
    print(round(classif,2))
    cat(" \n")
    
    cat("Posterior classification table:", "\n")
    cat("     --> mean of posterior probabilities in each class", "\n")
    print(round(cl.table,4))

cl.3 <- cl.table
entropy <- sum(mix.prob*log(mix.prob))
(ICL1_k3 <- BIC3 - entropy)
(ICL2_k3 <- BIC3 -2*sum(log((apply(mix.prob, MARGIN = 1, FUN = max)))))

entropy_3 <- 1 - (-sum(rowSums(mix.prob*log(mix.prob)))/(nrow(result$data_wide)*log(ncol(mix.prob))))
stat_i <- matrix(ncol=3, nrow=398)
cumHaz <- -log(result2$group_spec_surv)
stat <- 0
  for(i in 1:nrow(result$data_wide)){
    sums <- 0
    for(j in 1:ncol(result$theta$pj)){
    # i <- j <- 1
    big <- result$theta$ranef[[j]][i,]
    sum <- mix.prob[i,j]*(data_wide$delta[i]-cumHaz[i,j])*big
    sums <- sum+sums
    # stat <- sum+stat
    }
  stat_i[i,] <- sums
  }
# mean_stat <- tcrossprod(stat)/n
n <- nrow(result$data_wide)
n <- nrow(result$data_wide)
var <- matrix(0, nrow=3, ncol=3)
for(i in 1:nrow(result$data_wide)){
    var <- var + (tcrossprod(stat_i[i,]))
}
var <- var-tcrossprod(colSums(stat_i))/n
(score_stat <- t(colSums(stat_i))%*%matrixcalc::matrix.inverse(var)%*%colSums(stat_i))
(scoretest_pvalue3 <- 1-pchisq(score_stat, 3))
# score_mean <- sum()
# result2$predicted_data |> group_by(cluster)
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = data_wide, mix.prob = mix.prob, survform = survform)
library(SurvMetrics)
it3 <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)

rownames(mix.prob) <- result$data_wide$ID
yhat <- pred_yis(data_long=result$data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
preds <- rowSums(yhat * mix.prob[as.character(result$data_long$ID),])
local_MSE_y_in_3 <- Metrics::mse(actual = result$data_long$y, predicted = preds)


#! testing predictions
  n <- nrow(test_data_wide)
  # m <- nrow(betas_matrix)
  x <- model.frame(gform, data_wide)
# mixing proportions using nnet

    # For more than 2 components, use multinomial
    # postprob <- e_step_result$postprob

    nn_model <- nnet::nnet(
      x = x,          # Input features
      y = mix.prob,          # Target is component probabilities
      size = 10,      # Hidden layer size
      decay = 0.01,   #! Weight decay: this is something which can be tweaked
      maxit = 5000,   # Max iterations
      trace = FALSE,  # No output during training
      softmax = TRUE,  # Use softmax for multinomial
      entropy=TRUE # maximum conditional likelihood fitting
    )
    
    # Predict probabilities (this code doesn't work for some stupid reason, 
  # extract directly from the nnet object)
x_new <-  model.frame(gform, test_data_wide)
    prob <- predict(nn_model, x_new)
    # new_prior <- fitted(nn_model)
estep.result <- estep(prior_mixing=prob, theta=result$theta, long_form_o, long_form_g, ranform, survform,data_long =  test_data_long,data_wide =  test_data_wide)
mix.prob <- estep.result$P
 result2 <- plot_survival_curves(
  data_wide = test_data_wide,
  mix_prob = mix.prob,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "Out-Sample NN-JLCM"
)
nnjm_test <- result2
print(result2$plot)
# mix.prob <- estep.result$P
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = test_data_wide, mix.prob = mix.prob, survform = survform)

it3_test <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)


rownames(mix.prob) <- test_data_wide$ID

ranef <- pred_ranef(data_long=test_data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
result$theta$ranef <- ranef
yhat <- pred_yis(data_long=test_data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
preds <- rowSums(yhat * mix.prob[as.character(test_data_long$ID),])
local_MSE_y_out_3 <- Metrics::mse(actual = test_data_long$y, predicted = preds)

load("NNEM/Application/nnjlcmk4.RData")
# load("NNEM/Application/sd_res_Td.RData")
# result$theta_sd$`Survival Submodel` <- sd_res_Td

k <- nrow(result$theta_sd$`Longitudinal Submodel`)+nrow(result$theta_sd$`Survival Submodel`)
ll4 <- result$loglik
(BIC4=-2*ll4+k*log(nrow(result$data_wide)))
(AIC4=-2*ll4+2*k)


shape_params <- result$theta$shape
scale_params <-  result$theta$scale
gammas_params <-  result$theta$gamma_matrix
#! training predictions
# result$theta$shape <- shape_params
# result$theta$scale <- scale_params
# result$theta$gamma_matrix <- gammas_params
estep.result <- estep(prior_mixing=result$theta$pj, result$theta, long_form_o, long_form_g, ranform, survform, result$data_long, result$data_wide)
mix.prob <- result$postprob
# mix.prob <- estep.result$P
entropy <- sum(mix.prob*log(mix.prob))
(ICL1_k4 <- BIC4 - entropy)
(ICL2_k4 <- BIC4 -2*sum(log((apply(mix.prob, MARGIN = 1, FUN = max)))))

# Assuming you have extracted parameters from your other model:
result2 <- plot_survival_curves(
  data_wide = data_wide,
  mix_prob = result$postprob,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "NN-JLCM"
)
print(result2$plot)
pprob <- cbind(result$C, mix.prob)
    classif <- NULL
    classifY <- NULL
    cl.table <- NULL
    thr.table <- NULL
threshold=c(0.7,0.8,0.9)
ng <- ncol(mix.prob)
    for (g in 1:ng) 
    {
    temp<- subset(pprob,pprob[,1]==g)
    # tempY<- subset(pprobY,pprobY[,1]==g)
    temp1<-apply(temp[,1+1:ng],2,mean)
    # temp1Y<-apply(tempY[,1+1:ng],2,mean)
    cl.table<-rbind(cl.table,temp1)
    classif<-cbind(classif,length(temp[,1]))
    # classifY<-cbind(classifY,length(tempY[,1]))
    if(!is.null(threshold)) thr.table <- cbind(thr.table,sapply(threshold,function(x) length(which(temp[,1+g]>x)))/length(temp[,1]))  
   }
   
   classif <- rbind(classif,100*classif/nrow(result$data_wide))
  #  classifY <- rbind(classifY,100*classifY/x$ns)
      if(!is.null(thr.table))
   {
    thr.table <- 100*thr.table  
    
    rownames(thr.table) <- paste("prob>",threshold,sep="") 
    colnames(thr.table) <- paste("class",1:ng,sep="") 
   }
    if(!is.null(thr.table))
    {
     cat("Posterior probabilities above a threshold (%):", "\n")
     print(round(thr.table,2))
     cat(" \n")    
    }
   rownames(cl.table)<-paste("class",1:ng,sep="")
   colnames(cl.table)<-paste("prob",1:ng,sep="")
   colnames(classif)<-paste("class",1:ng,sep="")
   rownames(classif)<-c("N","%")
  #  colnames(classifY)<-paste("class",1:ng,sep="")
  #  rownames(classifY)<-c("N","%")
    cat(" \n")
    cat("Posterior classification based on longitudinal and time-to-event data:", "\n")
    print(round(classif,2))
    cat(" \n")
    
    cat("Posterior classification table:", "\n")
    cat("     --> mean of posterior probabilities in each class", "\n")
    print(round(cl.table,4))
entropy_4 <- 1 - (-sum(rowSums(mix.prob*log(mix.prob)))/(nrow(result$data_wide)*log(ncol(mix.prob))))
stat_i <- matrix(ncol=3, nrow=398)
cumHaz <- -log(result2$group_spec_surv)
stat <- 0
  for(i in 1:nrow(result$data_wide)){
    sums <- 0
    for(j in 1:ncol(result$theta$pj)){
    # i <- j <- 1
    big <- result$theta$ranef[[j]][i,]
    sum <- mix.prob[i,j]*(data_wide$delta[i]-cumHaz[i,j])*big
    sums <- sum+sums
    # stat <- sum+stat
    }
  stat_i[i,] <- sums
  }
# mean_stat <- tcrossprod(stat)/n
n <- nrow(result$data_wide)
var <- matrix(0, nrow=3, ncol=3)
for(i in 1:nrow(result$data_wide)){
    var <- var + (tcrossprod(stat_i[i,]))
}
var <- var-tcrossprod(colSums(stat_i))/n
(score_stat4 <- t(colSums(stat_i))%*%matrixcalc::matrix.inverse(var)%*%colSums(stat_i))
(scoretest_pvalue4 <- 1-pchisq(score_stat4, 3))


entropy_4 <- 1 - (-sum(rowSums(mix.prob*log(mix.prob)))/(nrow(result$data_wide)*log(ncol(mix.prob))))

# m_step_result <- mstep(postprob=mix.prob,theta, data_wide, data_long, ranform, long_form_g, long_form_o, gform, survform)
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = data_wide, mix.prob = mix.prob, survform = survform)
library(SurvMetrics)
it4 <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)

rownames(mix.prob) <- result$data_wide$ID
yhat <- pred_yis(data_long=result$data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
preds <- rowSums(yhat * mix.prob[as.character(result$data_long$ID),])
local_MSE_y_in_4 <- Metrics::mse(actual = result$data_long$y, predicted = preds)


#! testing predictions
  n <- nrow(test_data_wide)
  # m <- nrow(betas_matrix)
  x <- model.frame(gform, data_wide)
# mixing proportions using nnet

    # For more than 2 components, use multinomial
    # postprob <- e_step_result$postprob

    nn_model <- nnet::nnet(
      x = x,          # Input features
      y = mix.prob,          # Target is component probabilities
      size = 10,      # Hidden layer size
      decay = 0.01,   #! Weight decay: this is something which can be tweaked
      maxit = 5000,   # Max iterations
      trace = FALSE,  # No output during training
      softmax = TRUE,  # Use softmax for multinomial
      entropy=TRUE # maximum conditional likelihood fitting
    )
# ?nnet::nnet    
length(nn_model$wts)
    # Predict probabilities (this code doesn't work for some stupid reason, 
  # extract directly from the nnet object)
x_new <-  model.frame(gform, test_data_wide)
    prob <- predict(nn_model, x_new)
    # new_prior <- fitted(nn_model)
estep.result <- estep(prior_mixing=prob, theta=result$theta, long_form_o, long_form_g, ranform, survform,data_long =  test_data_long,data_wide =  test_data_wide)
mix.prob <- estep.result$P

 result2 <- plot_survival_curves(
  data_wide = test_data_wide,
  mix_prob = mix.prob,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "JLCM"
)
print(result2$plot)
# mix.prob <- estep.result$P
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = test_data_wide, mix.prob = mix.prob, survform = survform)

it4_test <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)


rownames(mix.prob) <- test_data_wide$ID

ranef <- pred_ranef(data_long=test_data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
result$theta$ranef <- ranef
yhat <- pred_yis(data_long=test_data_long, theta =  result$theta, long_form_g, long_form_o, ranform)
preds <- rowSums(yhat * mix.prob[as.character(test_data_long$ID),])
local_MSE_y_out_4 <- Metrics::mse(actual = test_data_long$y, predicted = preds)


df_aic <- data.frame(ng=c(2,3,4),loglik=c(ll1,ll3,ll4),BIC=c(BIC2,BIC3,BIC4),AIC=c(AIC2,AIC3,AIC4),measure_entropy=c(measure_entropy_2, entropy_3, entropy_4), ICL2=c(ICL2_k2, ICL2_k3, ICL2_k4), avg_diag=c(mean(diag(cl.2)), mean(diag(cl.3)), mean(diag(cl.table))), scoretest=c(scoretest_pvalue2, scoretest_pvalue3, scoretest_pvalue4), IBS_in=c(it2$IBS,it3$IBS,it4$IBS), MSE_in=c(local_MSE_y_in, local_MSE_y_in_3, local_MSE_y_in_4), IBS_out=c(it2_test$IBS, it3_test$IBS, it4_test$IBS), MSE_out=c(local_MSE_y_out_2, local_MSE_y_out_3, local_MSE_y_out_4))
df_aic 
# ?lcmm::summarytable
# Combine predictions with test data
dt <- cbind(yhat=yhat, test_data_long, mixprob=mix.prob[as.character(test_data_long$ID),])

# Calculate deciles
deciles_age65 <- quantile(result$data_long$age65, seq(0.1,1,0.1))
decile_breaks <- c(0, 0.566978, 0.813580, 1.055464, 1.248096, 1.442200,
                   1.648270, 1.851240, 2.076226, 2.366420, 3.308077)

# Rename columns for 3 clusters
setnames(dt, old = c("yhat.V1", "yhat.V2", "yhat.V3"),
             new = c("Yhat1", "Yhat2", "Yhat3"))
setnames(dt, old = c("mixprob.V1", "mixprob.V2", "mixprob.V3"),
             new = c("mixprob1", "mixprob2", "mixprob3"))

# Calculate weighted values for 3 clusters
dt[, `:=`(
  weighted_Y1 = mixprob1 * Yhat1,
  weighted_Y2 = mixprob2 * Yhat2,
  weighted_Y3 = mixprob3 * Yhat3
)]

# Create decile labels
decile_labels <- c("D1 (0-10%)", "D2 (10-20%)", "D3 (20-30%)", "D4 (30-40%)",
                   "D5 (40-50%)", "D6 (50-60%)", "D7 (60-70%)", "D8 (70-80%)",
                   "D9 (80-90%)", "D10 (90-100%)")

# Assign decile groups
dt[, decile_group := cut(measure_time,
                        breaks = decile_breaks,
                        labels = decile_labels,
                        include.lowest = TRUE)]

# Calculate means by decile for 3 clusters
decile_mean <- dt[, .(
  sum_weighted_Y1 = mean(weighted_Y1, na.rm = TRUE),
  sum_weighted_Y2 = mean(weighted_Y2, na.rm = TRUE),
  sum_weighted_Y3 = mean(weighted_Y3, na.rm = TRUE),
  count = .N
), by = decile_group]

print("Means by decile:")
print(decile_mean)

decile_mean[, `:=`(deciles = as.matrix(decile_breaks[-1]))]

# Save results
results_filename_nnem <- paste0("application_nnem_results.RData")
save(it2, it3, it4, file = results_filename_nnem)

# Plot 1: By decile number
plot3 <- decile_mean %>%
  pivot_longer(cols = starts_with("sum_weighted_Y"),
               names_to = "cluster",
               values_to = "sum_value") %>%
  mutate(cluster = gsub("sum_weighted_Y", "Cluster ", cluster),
         decile_numeric = as.numeric(gsub("D([0-9]+).*", "\\1", decile_group))) %>%
  ggplot(aes(x = decile_numeric, y = sum_value, color = cluster, group = cluster)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(title = "Trend of Weighted Y Values Across Deciles",
       x = "Decile",
       y = "Mean of Weighted Values",
       color = "Cluster") +
  scale_x_continuous(breaks = 1:10,
                     labels = unique(decile_mean$decile_group)[order(as.numeric(gsub("D([0-9]+).*", "\\1", unique(decile_mean$decile_group))))]) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot3)

# Show theta standard deviations
print(rbind(result$theta_sd$`Longitudinal Submodel`, result$theta_sd$`Survival Submodel`))

# Plot 2: By decile break values
plot3 <- decile_mean %>%
  pivot_longer(cols = starts_with("sum_weighted_Y"),
               names_to = "cluster",
               values_to = "sum_value") %>%
  mutate(cluster = gsub("sum_weighted_Y", "Cluster ", cluster),
         decile_numeric = as.numeric(gsub("D([0-9]+).*", "\\1", decile_group))) %>%
  ggplot(aes(x = deciles, y = sum_value, color = cluster, group = cluster)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(title = "Trend of Weighted Y Values Across Decile Breaks",
       x = "Decile Break Values",
       y = "Mean of Weighted Values",
       color = "Cluster") +
  scale_x_continuous(breaks = 0:10,
                     labels = round(t(decile_breaks), 3)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot3)

# Usage examples:

# # Example 1: Using JLCMM model
# load("NNEM/Application/init3model_application.RData")
# result1 <- plot_survival_curves(model = mj4b, data_wide = data_wide, title = "JLCMM", survform = survform)
# print(result1$plot)

# # Example 2: Using custom parameters (e.g., for NN-JLCM)
# # Assuming you have extracted parameters from your other model:
# result2 <- plot_survival_curves(
#   data_wide = data_wide,
#   mix_prob = result$postprob,
#   shape_params =shape_params,
#   scale_params = scale_params, 
#   gammas_params = gammas_params,
#   survform = survform,  # optional, will use default if not provided
#   title = "NN-JLCM"
# )
# print(result2$plot)

# Example 3: Compare both models side by side
# library(patchwork)
# result1$plot + result2$plot

#! JLCMM
#* in-sample
load("NNEM/Application/init2model_application.RData")
summary(mj2b)
m <- 2
# library(lme4)
# lmer(normMMSE ~ age65 + I(age65^2) + CEP + (age65|ID),paquidS)
sample_init <- estimates(mj2b, cholesky=F)

idWeib_shape <- which(grepl('Weibull2', names(sample_init)))
idWeib_scale <- which(grepl('Weibull1', names(sample_init)))

(scale_params <- exp(sample_init[idWeib_scale])^(-1/exp(sample_init[idWeib_shape])))
(shape_params <- exp(sample_init[idWeib_shape]))

# (scale_params <- 1/(sample_init[idWeib_scale]^(2)))
# (shape_params <- (sample_init[idWeib_shape])^2)
gammas_params <- sample_init[(max(idWeib_shape)+1):(max(idWeib_shape)+2)]

entropy <- 1 - (-sum(rowSums(mix.prob*log(mix.prob)))/(nrow(result$data_wide)*log(ncol(mix.prob))))

mix.prob <- as.matrix(mj2b$pprob[,c(3:4)])

ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = data_wide, mix.prob = mix.prob, survform = survform)
it_jm1 <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)
# v
#* out-sample
 result2 <- plot_survival_curves(
  data_wide = data_wide,
  mix_prob = mix.prob,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "JLCM"
)
jlcm_train <- result2
print(result2$plot)

mod <- mj2b
predy <- mod$pred$pred_ss
local_MSE_y_in_lcmm <- Metrics::mse(result$data_long$y, predy)
        # mod <- jlcmm
        est <- sample_init
          obj <- mod
          nclasses <- ncol(obj$pprob) - 2
          coefs <- obj$best
          coefend <- min(which(grepl('Weibull', names(coefs)))) - 1
          coefs_multilogit <- matrix(coefs[1:coefend], nrow = nclasses - 1)

          tmpX <- cbind(1, x)
          linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
          exps <- exp(cbind(linearval, 0))
          pi_logit <- exps / rowSums(exps)
          

          # Extract JLCMM parameters
          idxs_slope <- which(grepl(c('intercept', "age65 ", "measure_time2")[2], names(est)))

          b1 <- est[(idxs_slope-nclasses)]
          b2 <- est[idxs_slope]
          b3 <- est[(idxs_slope+nclasses)]
          df <- df2 <- data.frame(int=b1,mt=b2,mt2=b3)
          # df[2,] <- df2[1,]
          # df[1,] <- df2[3,]
          # df[3,] <- df2[4,]
          # df[4,] <- df[2,]
          local_sigma_est_lcmm <- est["stderr"]
          local_fixef_est_lcmm <- c(t(df))  # Based on your bias script
          local_shape_est_lcmm <- shape_params
          local_scale_est_lcmm <- scale_params
          
          # Variance parameters
          idxs_varcov <- which(grepl('varcov', names(est)))
          # idxs_varprop <- which(grepl('varprop', names(est)))
          local_chol_est_lcmm <- est[idxs_varcov]
          local_D_est_lcmm <- coefs[idxs_varcov]
          # local_omega_est_lcmm <- est[idxs_varprop]
          blank <- matrix(0, 3,3)
          diag(blank) <- local_D_est_lcmm[c(1,3,6)]
            blank[2,1] <- blank[1,2] <- local_D_est_lcmm[2]
            blank[3,1] <- blank[1,3] <- local_D_est_lcmm[4]
            blank[3,2] <- blank[2,3] <- local_D_est_lcmm[5]

              tmpX <- cbind(1, x_new)
              linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
              exps <- exp(cbind(linearval, 0))
              pi_logit_test <- exps / rowSums(exps)
              
              lcmm_theta <- result$theta
              lcmm_theta$beta_o <- est[min(idxs_varcov)-1]
              lcmm_theta$gamma_matrix <- gammas_params
              lcmm_theta$betas <- as.matrix(df)
              lcmm_theta$shape <- as.matrix(local_shape_est_lcmm)
              lcmm_theta$scale <- as.matrix(local_scale_est_lcmm)
              lcmm_theta$pj <- pi_logit
              lcmm_theta$D <- blank
              lcmm_theta$cholesky <- t(chol(blank)) 

              pred_mixing <- estep(prior_mixing=pi_logit_test, theta=lcmm_theta, long_form_o, long_form_g, ranform, survform,data_long =  test_data_long,data_wide =  test_data_wide)
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = test_data_wide, mix.prob = pred_mixing$P, survform = survform)

it2_test_lcmm <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)

summarytable(mj2b, which = c("entropy", "scoretest"))

rownames(pred_mixing$P) <- test_data_wide$ID
              predy_test_raw <- rowSums(pred_mixing$P[as.character(test_data_long$ID),] * 
                                        lcmm::predictY(obj, newdata=test_data_long)$pred)
              local_MSE_y_out_lcmm_2 <- Metrics::mse(test_data_long$y, predy_test_raw)

pprob2 <- postprob(mj2b)
prob2 <- mean(diag(pprob2[[2]]))

load("NNEM/Application/init3model_application.RData")
# summary(mj4b)
m <- 3
# library(lme4)
# lmer(normMMSE ~ age65 + I(age65^2) + CEP + (age65|ID),paquidS)
sample_init <- estimates(mj3b, cholesky=F)

idWeib_shape <- which(grepl('Weibull2', names(sample_init)))
idWeib_scale <- which(grepl('Weibull1', names(sample_init)))


(scale_params <- exp(sample_init[idWeib_scale])^(-1/exp(sample_init[idWeib_shape])))
(shape_params <- exp(sample_init[idWeib_shape]))


# (scale_params <- 1/(sample_init[idWeib_scale]^(2)))
# (shape_params <- (sample_init[idWeib_shape])^2)
gammas_params <- sample_init[(max(idWeib_shape)+1):(max(idWeib_shape)+2)]

mix.prob <- as.matrix(mj3b$pprob[,-c(1:2)])
 result2 <- plot_survival_curves(
  data_wide = data_wide,
  mix_prob = mix.prob,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "In-Sample JLCM"
)

jlcm_train <- result2
print(result2$plot)
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = data_wide, mix.prob = mix.prob, survform = survform)
it_jm3 <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)


mod <- mj3b
predy <- mod$pred$pred_ss
local_MSE_y_in_lcmm_3 <- Metrics::mse(result$data_long$y, predy)

        # mod <- jlcmm
        est <- sample_init
          obj <- mod
          nclasses <- ncol(obj$pprob) - 2
          coefs <- obj$best
          coefend <- min(which(grepl('Weibull', names(coefs)))) - 1
          coefs_multilogit <- matrix(coefs[1:coefend], nrow = nclasses - 1)

          tmpX <- cbind(1, x)
          linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
          exps <- exp(cbind(linearval, 0))
          pi_logit <- exps / rowSums(exps)
          

          # Extract JLCMM parameters
          idxs_slope <- which(grepl(c('intercept', "age65 ", "measure_time2")[2], names(est)))

          b1 <- est[(idxs_slope-nclasses)]
          b2 <- est[idxs_slope]
          b3 <- est[(idxs_slope+nclasses)]
          df <- df2 <- data.frame(int=b1,mt=b2,mt2=b3)
          # df[2,] <- df2[1,]
          # df[1,] <- df2[3,]
          # df[3,] <- df2[4,]
          # df[4,] <- df[2,]
          local_sigma_est_lcmm <- est["stderr"]
          local_fixef_est_lcmm <- c(t(df))  # Based on your bias script
          local_shape_est_lcmm <- shape_params
          local_scale_est_lcmm <- scale_params
          
          # Variance parameters
          idxs_varcov <- which(grepl('varcov', names(est)))
          # idxs_varprop <- which(grepl('varprop', names(est)))
          local_chol_est_lcmm <- est[idxs_varcov]
          local_D_est_lcmm <- coefs[idxs_varcov]
          # local_omega_est_lcmm <- est[idxs_varprop]
          blank <- matrix(0, 3,3)
          diag(blank) <- local_D_est_lcmm[c(1,3,6)]
            blank[2,1] <- blank[1,2] <- local_D_est_lcmm[2]
            blank[3,1] <- blank[1,3] <- local_D_est_lcmm[4]
            blank[3,2] <- blank[2,3] <- local_D_est_lcmm[5]

              tmpX <- cbind(1, x_new)
              linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
              exps <- exp(cbind(linearval, 0))
              pi_logit_test <- exps / rowSums(exps)
              
              lcmm_theta <- result$theta
              lcmm_theta$beta_o <- est[min(idxs_varcov)-1]
              lcmm_theta$gamma_matrix <- gammas_params
              lcmm_theta$betas <- as.matrix(df)
              lcmm_theta$shape <- as.matrix(local_shape_est_lcmm)
              lcmm_theta$scale <- as.matrix(local_scale_est_lcmm)
              lcmm_theta$pj <- pi_logit
              lcmm_theta$D <- blank
              lcmm_theta$cholesky <- t(chol(blank)) 

              pred_mixing <- estep(prior_mixing=pi_logit_test, theta=lcmm_theta, long_form_o, long_form_g, ranform, survform,data_long =  test_data_long,data_wide =  test_data_wide)
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = test_data_wide, mix.prob = pred_mixing$P, survform = survform)

it3_test_lcmm <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)

summarytable(mj3b, which = c("entropy", "scoretest"))

 result2 <- plot_survival_curves(
  data_wide = test_data_wide,
  mix_prob = pred_mixing$P,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "Out-Sample JLCM"
)

jlcm_test <- result2
print(result2$plot)
rownames(pred_mixing$P) <- test_data_wide$ID
              predy_test_raw <- rowSums(pred_mixing$P[as.character(test_data_long$ID),] * 
                                        lcmm::predictY(obj, newdata=test_data_long)$pred)
              local_MSE_y_out_lcmm_3 <- Metrics::mse(test_data_long$y, predy_test_raw)

pprob3 <- postprob(mj3b)
prob3 <- mean(diag(pprob3[[2]]))


load("NNEM/Application/init4model_application.RData")
# summary(mj4b)
m <- 4
# library(lme4)
# lmer(normMMSE ~ age65 + I(age65^2) + CEP + (age65|ID),paquidS)
sample_init <- estimates(mj4b, cholesky=F)

idWeib_shape <- which(grepl('Weibull2', names(sample_init)))
idWeib_scale <- which(grepl('Weibull1', names(sample_init)))



(scale_params <- exp(sample_init[idWeib_scale])^(-1/exp(sample_init[idWeib_shape])))
(shape_params <- exp(sample_init[idWeib_shape]))


# (scale_params <- 1/(sample_init[idWeib_scale]^(2)))
# (shape_params <- (sample_init[idWeib_shape])^2)
gammas_params <- sample_init[(max(idWeib_shape)+1):(max(idWeib_shape)+2)]

mix.prob <- as.matrix(mj4b$pprob[,-c(1:2)])
 result2 <- plot_survival_curves(
  data_wide = data_wide,
  mix_prob = mix.prob,
  shape_params =shape_params,
  scale_params = scale_params, 
  gammas_params = gammas_params,
  survform = survform,  # optional, will use default if not provided
  title = "Out-Sample"
)
result2$plot
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = data_wide, mix.prob = mix.prob, survform = survform)
it_jm4 <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)


mod <- mj4b
predy <- mod$pred$pred_ss
local_MSE_y_in_lcmm_4 <- Metrics::mse(result$data_long$y, predy)

        # mod <- jlcmm
        est <- sample_init
          obj <- mod
          nclasses <- ncol(obj$pprob) - 2
          coefs <- obj$best
          coefend <- min(which(grepl('Weibull', names(coefs)))) - 1
          coefs_multilogit <- matrix(coefs[1:coefend], nrow = nclasses - 1)

          tmpX <- cbind(1, x)
          linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
          exps <- exp(cbind(linearval, 0))
          pi_logit <- exps / rowSums(exps)
          

          # Extract JLCMM parameters
          idxs_slope <- which(grepl(c('intercept', "age65 ", "measure_time2")[2], names(est)))

          b1 <- est[(idxs_slope-nclasses)]
          b2 <- est[idxs_slope]
          b3 <- est[(idxs_slope+nclasses)]
          df <- df2 <- data.frame(int=b1,mt=b2,mt2=b3)
          # df[2,] <- df2[1,]
          # df[1,] <- df2[3,]
          # df[3,] <- df2[4,]
          # df[4,] <- df[2,]
          local_sigma_est_lcmm <- est["stderr"]
          local_fixef_est_lcmm <- c(t(df))  # Based on your bias script
          local_shape_est_lcmm <- shape_params
          local_scale_est_lcmm <- scale_params
          
          # Variance parameters
          idxs_varcov <- which(grepl('varcov', names(est)))
          # idxs_varprop <- which(grepl('varprop', names(est)))
          local_chol_est_lcmm <- est[idxs_varcov]
          local_D_est_lcmm <- coefs[idxs_varcov]
          # local_omega_est_lcmm <- est[idxs_varprop]
          blank <- matrix(0, 3,3)
          diag(blank) <- local_D_est_lcmm[c(1,3,6)]
            blank[2,1] <- blank[1,2] <- local_D_est_lcmm[2]
            blank[3,1] <- blank[1,3] <- local_D_est_lcmm[4]
            blank[3,2] <- blank[2,3] <- local_D_est_lcmm[5]

              tmpX <- cbind(1, x_new)
              linearval <- as.matrix(tmpX) %*% t(coefs_multilogit)
              exps <- exp(cbind(linearval, 0))
              pi_logit_test <- exps / rowSums(exps)
              
              lcmm_theta <- theta
              lcmm_theta$beta_o <- est[min(idxs_varcov)-1]
              lcmm_theta$gamma_matrix <- gammas_params
              lcmm_theta$betas <- as.matrix(df)
              lcmm_theta$shape <- as.matrix(local_shape_est_lcmm)
              lcmm_theta$scale <- as.matrix(local_scale_est_lcmm)
              lcmm_theta$pj <- pi_logit
              lcmm_theta$D <- blank
              lcmm_theta$cholesky <- t(chol(blank)) 

              pred_mixing <- estep(prior_mixing=pi_logit_test, theta=lcmm_theta, long_form_o, long_form_g, ranform, survform,data_long =  test_data_long,data_wide =  test_data_wide)
ret <- get_pred_surv(shape_params, scale_params, gammas_params = gammas_params, test_data = test_data_wide, mix.prob = pred_mixing$P, survform = survform)

it4_test_lcmm <- IBS(object=ret$brier_surv, sp_matrix=ret$Shat, IBSrange=ret$ti, n_cores=19,plot=TRUE, return_brier = TRUE)



rownames(pred_mixing$P) <- test_data_wide$ID
              predy_test_raw <- rowSums(pred_mixing$P[as.character(test_data_long$ID),] * 
                                        lcmm::predictY(obj, newdata=test_data_long)$pred)
              local_MSE_y_out_lcmm_4 <- Metrics::mse(test_data_long$y, predy_test_raw)

pprob4 <- postprob(mj4b)
prob4 <- mean(diag(pprob4[[2]]))

tmp <- data.frame(avg_diag=c(prob2, prob3, prob4),IBS_in=c(it_jm1$IBS, it_jm3$IBS, it_jm4$IBS), MSE_in=c(local_MSE_y_in_lcmm, local_MSE_y_in_lcmm_3, local_MSE_y_in_lcmm_4),IBS_out=c(it2_test_lcmm$IBS, it3_test_lcmm$IBS, it4_test_lcmm$IBS), MSE_out=c(local_MSE_y_out_lcmm_2, local_MSE_y_out_lcmm_3, local_MSE_y_out_lcmm_4))
df1 <- summarytable(mj2b,mj3b,mj4b,which=c("npm","loglik","BIC","AIC", "entropy", "scoretest", "ICL2"))
new <- cbind(df1,tmp)
new
df_aic
df
new
df2 <- summaryplot(mj1,mj2b,mj3b,mj4b,which=c("AIC","BIC","entropy"))
results_filename_jlcmm <- paste0("application_jlcmm_results.RData")
save(it_jm1,it_jm2, it_jm4, file = results_filename_jlcmm)


load("IBS_srem_test.RData"); load("IBS_srem.RData")
load("inout_jlct_mse.RData")
load("JLCT_ibs_test.RData")
jlct_ibs1_test <- jlct_ibs1
load("JLCT_ibs_train.RData")
load("predSurvJM.RData")
yhat <- unlist(preds$fitted.y, use.names = F)
y <- unlist(preds$y)

yhat_test <- NULL;
for(i in 1:398){
  n_ytmp <- length(preds$y[[i]])
  yhat_test_tmp <- preds$fitted.y[[i]][1:n_ytmp]
  yhat_test <- c(yhat_test,yhat_test_tmp)}
sremMse <- Metrics::mse(train_data_long$normMMSE,yhat_test)


p1 <- compare_brier_scores_gg(model1_results = it3,model2_results = it_jm3)
p2 <- compare_brier_scores_gg(model1_results = it3_test,model2_results = it3_test_lcmm, inout="Out-Sample")
library(patchwork)
combined_plot_patchwork <- p1/p2 +
  plot_annotation(
    title = "",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )
combined_plot_patchwork
print(p2)

# With custom styling
p2 <- compare_brier_scores_gg(it4, it_jm4,
                             colors = c("#1f77b4", "#ff7f0e"))
print(p2)


library(survival)
par(mfrow=c(1,1))


jlcm_train$plot / nnjm_train$plot 
jlcm_test$plot / nnjm_test$plot

result$postprob <- estep.result$P
probYT1 <- mj2b$pprob$probYT1

event_time <- result$data_wide$event_time
delta <- result$data_wide$delta
agedem <- paquidS$agedem
age_init <- paquidS$age_init
delta <- paquidS$dem

cols <- c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#592E83", 
                "#048A81", "#54632B", "#B5532B","#1f77b4")
plot(mj4b, which="survival", bty="l")
lines(survfit(Surv(event_time,delta)~1, weights=mj4b$pprob$probYT1, data=data_wide), bty="l")
lines(survfit(Surv(event_time,delta)~1, weights=mj4b$pprob$probYT2, data=data_wide), bty="l", col=cols[4])
lines(survfit(Surv(event_time,delta)~1, weights=mj4b$pprob$probYT3, data=data_wide), bty="l",  col=cols[6])
lines(survfit(Surv(event_time,delta)~1, weights=mj4b$pprob$probYT4, data=data_wide), bty="l", cols=cols[1])


brier

test_data=data_wide

est <- estimates(mj2b)
# 
# est <- 
D <- sample_init[1,which(grepl('varcov', colnames(est)))]
cholesky <- metafor::vec2mat(D, diag=T)
cholesky[3,1] <-cholesky[2,2] 
cholesky[2,2] <- cholesky[1,3]
cholesky[1,3] <- cholesky[3,1]
# cholesky <- round(cholesky)
# init[1:2] <- exp(init[1:2])

scores <- numeric(n)
for(i in 1:n){
    scores[i] <- sbrier(brier_surv, Shat[i,], btime=event_time[i])
}


# alphas

# betas_matrix <- 
# betas_matrix <- round(matrix(init, ncol = 3))
idWeib_shape <- which(grepl('Weibull2', names(est)))
idWeib_scale <- which(grepl('Weibull1', names(est)))

(scale_params <- 1/(est[idWeib_scale])^2)
(shape_params <- (est[idWeib_shape])^2)
# scales_est <- init[c(10,8,6,4)]^(-1/init[c(11,9,7,5)]) 

# scales_est <- (init[,1]^2)^(-1/(init[,2]^2))
# shapes_est <- exp(sample_init[,2])

gammas_matrix <- est[(max(idWeib_shape)+1):(max(idWeib_shape)+2)]

# betas_matrix <- sample_init[,2:8]

mix.prob <- as.matrix(mj2b$pprob[,-c(1,2)])


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


p1 <- compare_brier_scores_gg(it_1, jlcmm_results)
print(p1)

# With custom styling
p2 <- compare_brier_scores_gg(it_1, jlcmm_results,
                             model1_name = "NN-JLCM", 
                             model2_name = "JLCMM",
                             colors = c("#1f77b4", "#ff7f0e"))
print(p2)

bsurv <- survfit(brier_surv~1)
ybar <- 1-summary(bsurv, time=ti[2])$surv
rfit <- coxph(Surv(event_time, delta) ~ CEP+male,data_wide)
psurv <- survfit(rfit, newdata= data_wide)
dim(psurv)
# data
# 2982

 yhat <- 1- Shat[,2]

# weights for brier score
 wt4 <- rttright(Surv(event_time, delta) ~ 1, times =ti[2], data_wide)
table(wt4 ==0)

brier1 <- sum(wt4 * (data_wide$delta - yhat)^2)/ sum(wt4)
brier0 <- sum(wt4 * (data_wide$delta - ybar)^2) / sum(wt4)
r2 <- 1- (brier1/brier0)
temp <- c(numerator= brier1, denominator = brier0, rsquared = r2)

library(SurvMetrics)


IBS(brier_surv, Shat[,median(ti)>ti], ti[median(ti)>ti])
# IBS(brier_surv,Shat)
# IAEISE(brier_surv, Shat[,median(ti)>ti], ti[median(ti)>ti])

# Cindex(brier_surv, Shat[,490], t_star = ti[490])

#' Compare Brier Scores Between Two Models (ggplot2 version)
#'
#' Creates a comparative plot of Brier scores over time for two models using ggplot2
#'
#' @param model1_results List containing IBS, brier_scores, and time_points for first model
#' @param model2_results List containing IBS, brier_scores, and time_points for second model  
#' @param model1_name Character string for first model name (default: "NN-JLCM")
#' @param model2_name Character string for second model name (default: "JLCMM")
#' @param colors Vector of two colors for the models (default: c("#2E86AB", "#A23B72"))
#' @param add_difference Logical, whether to add a subplot showing the difference (default: FALSE)
#' @param save_plot Logical, whether to save the plot (default: FALSE)
#' @param filename Character string for saved plot filename (default: "brier_comparison.png")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 6)
#' @param dpi Resolution for saved plots (default: 300)
#'
#' @return Returns the ggplot object for further customization
#'
#' @examples
#' # Basic comparison
#' p <- compare_brier_scores_gg(it_1, jlcmm_results)
#' 
#' # With custom names and colors
#' p <- compare_brier_scores_gg(it_1, jlcmm_results, 
#'                             model1_name = "My NN-JLCM", 
#'                             model2_name = "Standard JLCMM",
#'                             colors = c("#1f77b4", "#ff7f0e"))
#'
#' # With difference plot
#' p <- compare_brier_scores_gg(it_1, jlcmm_results, add_difference = TRUE)
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export
compare_brier_scores_gg <- function(model1_results, model2_results, 
                                   model1_name = "NN-JLCM", model2_name = "JLCMM", inout="In-Sample",
                                   colors = c("#2E86AB", "#A23B72"), add_difference = FALSE,
                                   save_plot = FALSE, filename = "brier_comparison.png",
                                   width = 10, height = 6, dpi = 300) {
  
  # Load required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }
  if (add_difference && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package is required for difference plots")
  }
  
  library(ggplot2)
  library(dplyr)
  if (add_difference) library(patchwork)
  
  # Extract data from results
  time1 <- model1_results$time_points
  brier1 <- model1_results$brier_scores
  ibs1 <- model1_results$IBS
  
  time2 <- model2_results$time_points  
  brier2 <- model2_results$brier_scores
  ibs2 <- model2_results$IBS
  
  # Create combined data frame
  df1 <- data.frame(
    time = time1,
    brier = brier1,
    model = model1_name,
    ibs = ibs1,
    stringsAsFactors = FALSE
  )
  
  df2 <- data.frame(
    time = time2,
    brier = brier2,
    model = model2_name,
    ibs = ibs2,
    stringsAsFactors = FALSE
  )
  
  combined_df <- rbind(df1, df2)
  
  # Create IBS reference data
  ibs_df <- data.frame(
    model = c(model1_name, model2_name),
    ibs = c(ibs1, ibs2),
    color = colors,
    stringsAsFactors = FALSE
  )
  
  # Main comparison plot
  p1 <- ggplot(combined_df, aes(x = time, y = brier, color = model)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_hline(data = ibs_df, aes(yintercept = ibs, color = model), 
               linetype = "dashed", size = 0.8, alpha = 0.7) +
    scale_color_manual(values = setNames(colors, c(model2_name, model1_name)),
                       labels = paste0(c(model2_name, model1_name), 
                                     " (IBS = ", round(c(ibs2, ibs1), 4), ")")) +
    labs(
      title = paste(inout,"Brier Score Comparison:", model1_name, "vs", model2_name),
      subtitle = paste("Lower values indicate better predictive performance"),
      x = "Time",
      y = "Brier Score",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", size = 0.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  # Add difference plot if requested
  if (add_difference) {
    # Determine common time range
    time_min <- max(min(time1, na.rm = TRUE), min(time2, na.rm = TRUE))
    time_max <- min(max(time1, na.rm = TRUE), max(time2, na.rm = TRUE))
    
    # Filter to common range
    df1_filtered <- df1 %>% filter(time >= time_min & time <= time_max)
    df2_filtered <- df2 %>% filter(time >= time_min & time <= time_max)
    
    # Interpolate to common time points if needed
    if (nrow(df1_filtered) != nrow(df2_filtered) || !all(df1_filtered$time == df2_filtered$time)) {
      common_times <- seq(time_min, time_max, length.out = min(nrow(df1_filtered), nrow(df2_filtered)))
      brier1_interp <- approx(df1_filtered$time, df1_filtered$brier, common_times)$y
      brier2_interp <- approx(df2_filtered$time, df2_filtered$brier, common_times)$y
    } else {
      common_times <- df1_filtered$time
      brier1_interp <- df1_filtered$brier
      brier2_interp <- df2_filtered$brier
    }
    
    difference_df <- data.frame(
      time = common_times,
      difference = brier1_interp - brier2_interp
    )
    
    mean_diff <- mean(difference_df$difference, na.rm = TRUE)
    
    # Difference plot
    p2 <- ggplot(difference_df, aes(x = time, y = difference)) +
      geom_line(color = "#2E8B57", size = 1.2) +
      geom_point(color = "#2E8B57", size = 1.2, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
      geom_hline(yintercept = mean_diff, linetype = "dotted", color = "#2E8B57", size = 0.8) +
      labs(
        title = paste("Brier Score Difference:", model1_name, "−", model2_name),
        subtitle = paste("Mean Difference =", round(mean_diff, 4), 
                        ifelse(mean_diff > 0, paste("(", model2_name, "better)"), 
                              paste("(", model1_name, "better)"))),
        x = "Time",
        y = paste("Difference (", model1_name, "−", model2_name, ")")
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", size = 0.5)
      )
    
    # Combine plots using patchwork
    final_plot <- p1 / p2 + plot_layout(heights = c(2, 1))
    
  } else {
    final_plot <- p1
  }
  
  # Save plot if requested
  if (save_plot) {
    ggsave(filename, final_plot, width = width, height = height, dpi = dpi, bg = "white")
    cat("Plot saved as:", filename, "\n")
  }
  
  return(final_plot)
}

#' Compare Multiple Models' Brier Scores (ggplot2 version)
#'
#' Creates a comparative plot for multiple models using ggplot2
#'
#' @param model_results_list List of model results, each containing IBS, brier_scores, time_points
#' @param model_names Vector of model names
#' @param colors Vector of colors for each model (optional)
#' @param save_plot Logical, whether to save the plot
#' @param filename Filename for saved plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution for saved plots
#' @param show_ibs_lines Logical, whether to show IBS reference lines (default: TRUE)
#' @param facet_by_model Logical, whether to create separate facets for each model (default: FALSE)
#'
#' @examples
#' # Compare three models
#' p <- compare_multiple_brier_scores_gg(
#'   list(it_1, jlcmm_results, other_model_results),
#'   c("NN-JLCM", "JLCMM", "Other Model")
#' )
#'
#' # With faceting
#' p <- compare_multiple_brier_scores_gg(
#'   list(it_1, jlcmm_results, other_model_results),
#'   c("NN-JLCM", "JLCMM", "Other Model"),
#'   facet_by_model = TRUE
#' )
#'
#' @import ggplot2
#' @import dplyr
#' @export
compare_multiple_brier_scores_gg <- function(model_results_list, model_names,
                                            colors = NULL, save_plot = FALSE,
                                            filename = "multiple_brier_comparison.png",
                                            width = 12, height = 7, dpi = 300,
                                            show_ibs_lines = TRUE, facet_by_model = FALSE, inout) {
  
  # Load required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }
  
  library(ggplot2)
  library(dplyr)
  
  n_models <- length(model_results_list)
  
  if (n_models > 8) {
    stop("Maximum 8 models supported for clear visualization")
  }
  
  if (is.null(colors)) {
    colors <- c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#592E83", 
                "#048A81", "#54632B", "#B5532B")[1:n_models]
  }
  
  if (length(model_names) != n_models) {
    stop("Number of model names must match number of models")
  }
  
  # Create combined data frame
  combined_df <- data.frame()
  ibs_values <- numeric(n_models)
  
  for (i in 1:n_models) {
    df_temp <- data.frame(
      time = model_results_list[[i]]$time_points,
      brier = model_results_list[[i]]$brier_scores,
      model = model_names[i],
      stringsAsFactors = FALSE
    )
    combined_df <- rbind(combined_df, df_temp)
    ibs_values[i] <- model_results_list[[i]]$IBS
  }
  
  # Create IBS reference data
  ibs_df <- data.frame(
    model = model_names,
    ibs = ibs_values,
    stringsAsFactors = FALSE
  )
  
  combined_df$model <- factor(combined_df$model, levels = model_names)
ibs_df$model <- factor(ibs_df$model, levels = model_names)
  # Create the plot
  p <- ggplot(combined_df, aes(x = time, y = brier, color = model)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 1.2, alpha = 0.6) +
    scale_color_manual(values = setNames(colors, model_names),
                       labels = paste0(model_names, " (IBS = ", round(ibs_values, 4), ")")) +
    labs(
      title = paste0(inout," Brier Score Comparison"),
      x = "Time",
      y = "Brier Score",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", size = 0.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2))
  
  # Add IBS reference lines if requested
  if (show_ibs_lines) {
    p <- p + geom_hline(data = ibs_df, aes(yintercept = ibs, color = model), 
                        linetype = "dashed", size = 0.6, alpha = 0.7)
  }
  
  # Add faceting if requested
  if (facet_by_model) {
    p <- p + facet_wrap(~model, scales = "free_y", ncol = ceiling(sqrt(n_models))) +
      theme(
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "none"  # Remove legend when faceting since it's redundant
      )
  }
  
  # Save plot if requested
  if (save_plot) {
    ggsave(filename, p, width = width, height = height, dpi = dpi, bg = "white")
    cat("Plot saved as:", filename, "\n")
  }
  
  return(p)
}

#' Create Summary Table of Model Performance
#'
#' Creates a formatted table comparing IBS values across models
#'
#' @param model_results_list List of model results
#' @param model_names Vector of model names
#' @param save_table Logical, whether to save as CSV
#' @param filename Filename for saved table
#'
#' @return Data frame with model comparison
#' @import dplyr
#' @export
create_model_summary_table <- function(model_results_list, model_names, 
                                      save_table = FALSE, filename = "model_summary.csv") {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }
  
  library(dplyr)
  
  # Extract IBS values and calculate statistics
  summary_df <- data.frame(
    Model = model_names,
    IBS = sapply(model_results_list, function(x) round(x$IBS, 4)),
    Min_Brier = sapply(model_results_list, function(x) round(min(x$brier_scores, na.rm = TRUE), 4)),
    Max_Brier = sapply(model_results_list, function(x) round(max(x$brier_scores, na.rm = TRUE), 4)),
    Mean_Brier = sapply(model_results_list, function(x) round(mean(x$brier_scores, na.rm = TRUE), 4)),
    SD_Brier = sapply(model_results_list, function(x) round(sd(x$brier_scores, na.rm = TRUE), 4)),
    stringsAsFactors = FALSE
  ) %>%
    arrange(IBS) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, Model, IBS, Mean_Brier, Min_Brier, Max_Brier, SD_Brier)
  
  if (save_table) {
    write.csv(summary_df, filename, row.names = FALSE)
    cat("Summary table saved as:", filename, "\n")
  }
  
  return(summary_df)
}


p1 <- compare_multiple_brier_scores_gg(list(it_jm3, it3, it_srem, jlct_ibs1), c('JLCM', 'NN-JLCM', 'SREM', 'JLCT'), inout="In-Sample")
p2 <- compare_multiple_brier_scores_gg(list(it3_test_lcmm, it3_test, it_srem_test, jlct_ibs1_test), c('JLCM', 'NN-JLCM', 'SREM', 'JLCT'), inout="Out-Sample")
p1/p2
