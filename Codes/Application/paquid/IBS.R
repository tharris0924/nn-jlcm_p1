#' IBS
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
#'
#' @return The integration of brierscore
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
#' # the default time points
#' IBS(object, sp_matrix, time_interest)
#' # a time range
#' IBS(object, sp_matrix, c(18:100))
#'
#' @importFrom survival Surv
#'
#'
#' @export
IBS <- function(object, sp_matrix, IBSrange = c(-2, -1)) {
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
    sp_matrix <- predictSurvProb.survreg(obj, test_data, distime)
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
    t_brier <- rep(0)
    for (i in 1:length(IBSrange)) {
      pre_sp <- sp_matrix[, i]
      t_star <- IBSrange[i]
      t_brier[i] <- Brier(object, pre_sp, t_star)
    }
    t_IBS <- 0
    for (i in 1:(length(IBSrange) - 1)) {
      t_IBS <- t_IBS + (IBSrange[i + 1] - IBSrange[i]) * t_brier[i]
    }
    t_IBS <- t_IBS / (range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) <- "IBS"
    return(round(t_IBS, 6))
  } else {
    t_brier <- rep(0)
    t_IBSrange <- range(IBSrange)
    p <- ncol(sp_matrix)
    # Simulation to generate discrete time series
    IBSrange <- seq(t_IBSrange[1], t_IBSrange[2], length = p)
    for (i in 1:length(IBSrange)) {
      pre_sp <- sp_matrix[, i]
      t_star <- IBSrange[i]
      t_brier[i] <- Brier(object, pre_sp, t_star)
    }
    t_IBS <- 0
    for (i in 1:(length(IBSrange) - 1)) {
      t_IBS <- t_IBS + (IBSrange[i + 1] - IBSrange[i]) * t_brier[i]
    }
    t_IBS <- t_IBS / (range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) <- "IBS"
    return(round(t_IBS, 6))
  }
}


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
#'
#' @return The integration of brierscore
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
#' # the default time points with parallel processing
#' IBS(object, sp_matrix, time_interest, parallel = TRUE)
#' # a time range with parallel processing
#' IBS(object, sp_matrix, c(18:100), parallel = TRUE, n_cores = 4)
#'
#' @importFrom survival Surv
#' @importFrom parallel makeCluster stopCluster parLapply detectCores clusterEvalQ clusterExport
#'
#'
#' @export
IBS <- function(object, sp_matrix, IBSrange = c(-2, -1), parallel = TRUE, n_cores = NULL) {
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
    sp_matrix <- predictSurvProb.survreg(obj, test_data, distime)
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

  # Helper function for parallel Brier score computation
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
    return(round(t_IBS, 6))
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
    return(round(t_IBS, 6))
  }
}


it <- IBS(object=brier_surv, sp_matrix=Shat, IBSrange=ti, n_cores=19,plot=TRUE, return_brier = TRUE)
it <- IBS(brier_surv, Shat[,median(ti)>ti], ti[median(ti)>ti], n_cores=19,plot=TRUE, return_brier = TRUE)
