# density of survival times with censoring
  f_td <-  function(shape, scale, t, d) {
    ((dweibull(t,shape=shape, scale = scale)^d)) * (pweibull(t, shape = shape, scale =scale,lower.tail = F)^((1-d))) 
  }


survform=~CEP+male
# gammas_params <- theta$gamma_matrix
density_weibull <- function(shape_params, scale_params, gammas_params, data_wide, survform) {

  survform <- update(survform, ~.-1)
  W <- model.matrix(survform, data_wide)
  N <- nrow(data_wide)
  k <- nrow(as.matrix(shape_params))


  f_tds <- matrix(0,nrow = N, ncol = k)
  # w <- model.matrix(~ -1+X3+X4+X5, data=data_wide)
  delta <- data_wide$delta
  # initial density function for Weibull
  
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
    ti <- data_wide$event_time
    for(i in 1:N){
    f_tds[i, j] <- f_td(shape, new_scales[i], ti[i], delta[i])
    }
    }
  return(f_tds)
}

# ftds <- density_weibull(shape_params = theta$shape, scale_params = theta$scale, gammas_params = gammas_params, data_wide, ~CEP+male)
long_form_g=y~1+ measure_time +(measure_time2)
long_form_o=~-1+ CEP

ranform=~1+ measure_time +(measure_time2)
# data_long <- datalong
marginal_density_yis <- function(data_long, theta, long_form_g, long_form_o,ranform) {
  # parms prep
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

  N <- length(unique(data_long$ID)) #;print(N) # number of unique subjects under study.
  m <- nrow(betas_matrix) #;print(m);print(betas_matrix) # number of unique groups for mixture model
  f_yi <- matrix(0,nrow = N, ncol = m)
  nis<-table(data_long$ID)
  IDs<-unique(data_long$ID)
  # length(nis)

  D_matrices <- tcrossprod(theta$cholesky);# print(D_matrices)
  for (k in 1:m) {
    sigma_g <- sigmas # error variances 
    D_g <- D_matrices # random effect variance-covariance matrices
    beta_g <- betas_matrix[k,] # fixed effects
  
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
      Xi_beta_k <- matrix(X_ig %*% beta_g + X_io %*% beta_o) # Fixed Effects Linear Predictor
   
      q <- ncol(Z_ig)
      Y_i <- subset$y # Longitudinal recordings

      #? calculating the distributional parameters for the group specific random effects
      # D_g_inv <- solve(D_g)
      # 3-by-2 dimensions
      V_ig <- as.matrix(Z_ig %*% (D_g) %*% t(Z_ig) + sigma_g^2 * diag(nrow(X_ig)))
      # simgas[i] <- list(Sigma_i)
      f_yi[i, k] <- mvnfast::dmvn(t(Y_i), Xi_beta_k, V_ig)
    }
    }
  return(f_yi)
}


# example usage
ysu <- marginal_density_yis(
  data_long = datalong, 
theta=theta,
long_form_g = long_form_g,
long_form_o = long_form_o,
ranform = ranform)

#' IBS with Parallel Processing
#'
#' IBS is  an integrated version of the Brier which is used to calculate the integration of the Brier Score.
#' The Brier Score is the mean square difference between the true classes and the predicted probabilities.
#' Basically, the IBS is an integrated
#' weighted squared distance between the estimated survival function and the empirical
#' survival function. The inverse probability censoring weighting(IPCW) is used to adjust for censoring.
#'
#' The percentage of censored observations increases in time, and this will surely affect the dispersion of the empirical Brier Score.
#' The question of how censoring in finite samples acts on the distribution of our measures of inaccuracy is an interesting subject.
#' Our recommendation is to choose t* in a way that censoring is not too heavy (for example, the median follow-up time).
#' We also prefer measures with integrated loss functions since they will reflect inaccuracy over an interval rather than just at one point in time.
#' In addition, the corresponding empirical measures are likely to have lower dispersion, because censored observations contribute their estimated event-free probabilities to the integrand until the censoring occurs.
#'
#' @param object object of class \code{Surv} in the testing set created by Surv function.
#' @param sp_matrix a matrix or data.frame of predicted values of survival probabilities for the testing set.
#' @param IBSrange a vector contains all discrete time points corresponding to the predicted probability in sp_matrix.
#' Or the scale you want to get the IBS; and if it is a single point the return value will be the Brier Score at the timepoint.
#' @param parallel logical, whether to use parallel processing (default: TRUE)
#' @param n_cores integer, number of cores to use for parallel processing. If NULL, uses all available cores minus 1.
#' @param plot logical, whether to create a plot of Brier scores over time (default: FALSE)
#' @param return_brier logical, whether to return individual Brier scores along with IBS (default: FALSE)
#'
#' @return If return_brier = FALSE, returns the integration of brierscore (IBS). 
#' If return_brier = TRUE, returns a list containing:
#' \itemize{
#'   \item IBS: the integrated Brier score
#'   \item brier_scores: vector of individual Brier scores at each time point
#'   \item time_points: corresponding time points
#' }
#'
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#' HooraMoradian, DenisLarocque, & FranoisBellavance. (2017). \\(l_1\\) splitting rules in survival forests. Lifetime Data Analysis, 23(4), 671–691.
#'
#' Graf, Erika, Schmoor, Claudia, Sauerbrei, & Willi, et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. Statist. Med., 18(1718), 2529-2545.
#'
#' Brier, G. W. . (1950). Verification of forecasts expressed in terms of probability. Monthly Weather Review, 78.
#'
#' Gneiting, T. , &  Raftery, A. E. . (2007). Strictly Proper Scoring Rules, Prediction, and Estimation.

#' @examples
#' library(survival)
#' library(SurvMetrics)
#' library(parallel)
#' set.seed(123)
#' N <- 100
#' mydata <- SDGM4(N, p = 20, c_step = -0.5)
#' index.train <- sample(1:N, 2 / 3 * N)
#' data.train <- mydata[index.train, ]
#' data.test <- mydata[-index.train, ]
#'
#' time_interest <- sort(data.train$time[data.train$status == 1])
#' sp_matrix <- matrix(sort(runif(nrow(data.test) * length(time_interest)),
#'   decreasing = TRUE
#' ), nrow = nrow(data.test))
#' object <- Surv(data.test$time, data.test$status)
#'
#' # the default time points with parallel processing and plotting
#' result <- IBS(object, sp_matrix, time_interest, parallel = TRUE, plot = TRUE, return_brier = TRUE)
#' # a time range with parallel processing and custom cores
#' IBS(object, sp_matrix, c(18:100), parallel = TRUE, n_cores = 4, plot = TRUE)
#'
#' @importFrom survival Surv
#' @importFrom parallel makeCluster stopCluster parLapply detectCores clusterEvalQ clusterExport
#' @importFrom graphics plot lines points legend title
#' @importFrom grDevices dev.new
#'
#'
#' @export
IBS <- function(object, sp_matrix, IBSrange = c(-2, -1), parallel = TRUE, n_cores = NULL, plot = FALSE, return_brier = FALSE) {
  # case1、coxph AND testing set
  if (inherits(object, "coxph")) {
    obj <- object
    test_data <- sp_matrix

    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))

    mat_coxph <-
      pec::predictSurvProb(obj, test_data, distime) # get the survival probability matrix
    object_coxph <- Surv(test_data$time, test_data$status)

    object <- object_coxph
    sp_matrix <- mat_coxph
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
    } # the fixed time point
  }


  # case2、RSF AND testing set
  if (inherits(object, c("rfsrc"))) {
    obj <- object
    test_data <- sp_matrix

    mat_rsf <-
      predict(obj, test_data)$survival # get the survival probability matrix
    object_rsf <- Surv(test_data$time, test_data$status)

    object <- object_rsf
    sp_matrix <- mat_rsf
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
    } # the fixed time point
  }


  # case3 survreg AND testing set
  if (inherits(object, c("survreg"))) {
    obj <- object
    test_data <- sp_matrix

    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))

    object <- Surv(test_data$time, test_data$status)
    sp_matrix <- pec::predictSurvProb.survreg(obj, test_data, distime)
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
    } # the fixed time point
  }

  if (!is.numeric(IBSrange)) {
    stop("The class of the IBSrange must be numeric! or use the default setting")
  }

  if (any(is.na(IBSrange))) {
    stop("Cannot calculate IBS in the interval containing NA")
  }

  #default time range for IBS()
  if(max(IBSrange) <= 0){
    IBSrange = range(object[, 1])
  }

  if (!inherits(object, "Surv")) {
    stop("object is not of class Surv")
  }

  if (missing(object)) {
    stop("The survival object of the testing set is missing")
  }

  if (missing(sp_matrix)) {
    stop("The prediction of the survival probability matrix is missing")
  }

  if (any(is.na(object))) {
    stop("The input vector cannot have NA")
  }

  if (any(is.na(sp_matrix))) {
    stop("The input probability matrix cannot have NA")
  }

  if (length(object) != nrow(sp_matrix)) {
    stop(
      "number of rows of the sp_matrix and the survival object have different lengths"
    )
  }

  if (any(IBSrange <= 0)) {
    stop("The integration interval must be positive")
  }

  if (any(diff(IBSrange) <= 0)) {
    stop("The integral interval value must increase")
  }

  # Helper function for plotting Brier scores
  plot_brier_scores <- function(time_points, brier_scores, ibs_value) {
    # Create the plot
    plot(time_points, brier_scores, 
         type = "l", 
         col = "blue", 
         lwd = 2,
         xlab = "Time", 
         ylab = "Brier Score",
         main = paste("Brier Scores Over Time\n(IBS =", round(ibs_value, 4), ")"),
         cex.main = 1.2,
         cex.lab = 1.1)
    
    # Add points
    points(time_points, brier_scores, col = "blue", pch = 19, cex = 0.8)
    
    # Add horizontal line for mean Brier score (IBS approximation)
    mean_brier <- mean(brier_scores)
    abline(h = mean_brier, col = "red", lty = 2, lwd = 2)
    
    # Add grid for better readability
    grid(col = "lightgray", lty = 3)
    
    # Add legend
    legend("topright", 
           legend = c("Brier Score", "Mean Brier Score"),
           col = c("blue", "red"),
           lty = c(1, 2),
           lwd = c(2, 2),
           pch = c(19, NA),
           cex = 0.9,
           bg = "white")
  }
  # Helper function to return results
  format_results <- function(ibs_value, time_points = NULL, brier_scores = NULL, return_brier = FALSE, plot = FALSE) {
    if (plot) {
      plot_brier_scores(time_points, brier_scores, ibs_value)
    }
    
    if (return_brier) {
      return(list(
        IBS = ibs_value,
        brier_scores = brier_scores,
        time_points = time_points
      ))
    } else {
      return(ibs_value)
    }
  }
  compute_brier_parallel <- function(indices, object, sp_matrix, IBSrange) {
    sapply(indices, function(i) {
      pre_sp <- sp_matrix[, i]
      t_star <- IBSrange[i]
      Brier(object, pre_sp, t_star)
    })
  }

  # Determine the input dimension
  if (ncol(data.frame("sp" = sp_matrix)) == 1 &
    length(IBSrange) == 1) {
    bs <- Brier(object, sp_matrix, IBSrange)
    names(bs) <- "Brier Score"
    return(round(bs, 6))
  } else if ((ncol(data.frame("sp" = sp_matrix)) - 1) * (length(IBSrange) -
    1) == 0) {
    stop("The length is illegal")
  } else if (ncol(sp_matrix) == length(IBSrange)) {
    IBSrange <- sort(IBSrange)
    
    # Parallel computation of Brier scores
    if (parallel && length(IBSrange) > 1) {
      # Set up parallel processing
      if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
      }
      
      # Only use parallel processing if we have multiple cores and multiple time points
      if (n_cores > 1 && length(IBSrange) > 2) {
        cl <- parallel::makeCluster(n_cores)
        
        # Export necessary objects and functions to cluster
        parallel::clusterExport(cl, c("object", "sp_matrix", "IBSrange", "Brier"), 
                               envir = environment())
        
        # Load required packages on each worker
        parallel::clusterEvalQ(cl, {
          library(survival)
          library(SurvMetrics)
        })
        
        # Split indices for parallel processing
        indices_list <- split(1:length(IBSrange), 
                             cut(1:length(IBSrange), n_cores, labels = FALSE))
        
        # Compute Brier scores in parallel
        t_brier_list <- parallel::parLapply(cl, indices_list, 
                                          compute_brier_parallel,
                                          object = object, 
                                          sp_matrix = sp_matrix,
                                          IBSrange = IBSrange)
        
        # Stop cluster
        parallel::stopCluster(cl)
        
        # Combine results
        t_brier <- unlist(t_brier_list)
      } else {
        # Fall back to sequential processing
        t_brier <- sapply(1:length(IBSrange), function(i) {
          pre_sp <- sp_matrix[, i]
          t_star <- IBSrange[i]
          Brier(object, pre_sp, t_star)
        })
      }
    } else {
      # Sequential computation
      t_brier <- sapply(1:length(IBSrange), function(i) {
        pre_sp <- sp_matrix[, i]
        t_star <- IBSrange[i]
        Brier(object, pre_sp, t_star)
      })
    }
    
    t_IBS <- 0
    for (i in 1:(length(IBSrange) - 1)) {
      t_IBS <- t_IBS + (IBSrange[i + 1] - IBSrange[i]) * t_brier[i]
    }
    t_IBS <- t_IBS / (range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) <- "IBS"
    ret <- format_results(t_IBS,time_points = IBSrange, brier_scores=t_brier, return_brier = return_brier, plot=plot)
    return(ret)
  } else {
    t_IBSrange <- range(IBSrange)
    p <- ncol(sp_matrix)
    # Simulation to generate discrete time series
    IBSrange <- seq(t_IBSrange[1], t_IBSrange[2], length = p)
    
    # Parallel computation of Brier scores
    if (parallel && p > 1) {
      # Set up parallel processing
      if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
      }
      
      # Only use parallel processing if we have multiple cores and multiple time points
      if (n_cores > 1 && p > 2) {
        cl <- parallel::makeCluster(n_cores)
        
        # Export necessary objects and functions to cluster
        parallel::clusterExport(cl, c("object", "sp_matrix", "IBSrange", "Brier"), 
                               envir = environment())
        
        # Load required packages on each worker
        parallel::clusterEvalQ(cl, {
          library(survival)
          library(SurvMetrics)
        })
        
        # Split indices for parallel processing
        indices_list <- split(1:p, cut(1:p, n_cores, labels = FALSE))
        
        # Compute Brier scores in parallel
        t_brier_list <- parallel::parLapply(cl, indices_list, 
                                          compute_brier_parallel,
                                          object = object, 
                                          sp_matrix = sp_matrix,
                                          IBSrange = IBSrange)
        
        # Stop cluster
        parallel::stopCluster(cl)
        
        # Combine results
        t_brier <- unlist(t_brier_list)
      } else {
        # Fall back to sequential processing
        t_brier <- sapply(1:p, function(i) {
          pre_sp <- sp_matrix[, i]
          t_star <- IBSrange[i]
          Brier(object, pre_sp, t_star)
        })
      }
    } else {
      # Sequential computation
      t_brier <- sapply(1:p, function(i) {
        pre_sp <- sp_matrix[, i]
        t_star <- IBSrange[i]
        Brier(object, pre_sp, t_star)
      })
    }
    
    t_IBS <- 0
    for (i in 1:(length(IBSrange) - 1)) {
      t_IBS <- t_IBS + (IBSrange[i + 1] - IBSrange[i]) * t_brier[i]
    }
    t_IBS <- t_IBS / (range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) <- "IBS"

    return(list(t_IBS, t_brier, IBSrange))
  }
}