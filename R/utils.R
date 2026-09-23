# helper functions for the drcmd function

#' @title find_missing_pattern
#' @description Find the missing pattern in the data
#' @param Y A vector or data frame containing outcome values
#' @param A A vector or data frame  containing treatment variable values
#' @param X A data frame containing covariate values
#' @param W A data frame containing proxy variable values
#' @param min_complete Number of complete cases below which to issue a warning
#'
#' @return A character string containing the missing pattern
#' @keywords internal
#' @examples
#' \dontrun{
#' n <- 200
#' X <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, 0.5)
#' Y <- rnorm(n)
#' R <- rbinom(n, 1, 0.7)
#' Y[R == 0] <- NA
#' W <- X[, 0]
#' result <- find_missing_pattern(Y, A, X, W)
#' result$Z  # variables without missingness
#' result$U  # variables with missingness
#' result$R  # complete case indicator
#' }
find_missing_pattern <- function(Y,A,X,W,min_complete=10L) {

  # Combine variables into a single data frame
  data <- cbind(X,W,
                data.frame(Y=Y,A=A))

  # find variables that are never missing
  na_counts <- colSums(is.na(data))
  never_missing <- colnames(data)[na_counts == 0]
  if(length(never_missing) == 0) {
    stop("Error: drcmd requires data to have at least one variable that is never missing")
  }

  # find variables that are sometimes missing
  sometimes_missing <- colnames(data)[na_counts > 0]

  # make variable that indicates complete cases and df of all complete data
  R <- as.numeric(rowSums(is.na(data)) == 0)
  if (all(R == 0)) {
    stop("Error: drcmd requires data to have at least one complete case")
  }
  if (sum(R) < min_complete || mean(R) < 0.01) {
    warning('Only ', sum(R), ' complete cases are available. Results may be unstable')
  }

  Z <- data[, never_missing, drop=FALSE]

  # If 'Y' is in Z rename to 'y' (superlearner doesn't allow covariates named Y)
  if('Y' %in% colnames(Z)) {
    colnames(Z)[colnames(Z) == 'Y'] <- 'y'
  }

  X[] <- lapply(X, function(x) {
    if (!anyNA(x)) return(x)
    replacement <- if (is.numeric(x)) 0 else x[which(!is.na(x))[1L]]
    x[is.na(x)] <- replacement
    x
  })
  Y[is.na(Y)] <- 0
  A[is.na(A)] <- 0

  return(list(Z=Z,R=R,X=X,Y=Y,A=A,
              U=sometimes_missing))

}

#' @title check_r_ind
#'
#' @description Check if the missingness indicator R is defined correctly
#'
#' @param data A data frame
#' @param Y A character string containing outcome variable name
#' @param A A character string containing treatment variable name
#' @param X A character vector containing covariate variable names
#' @param W A character vector containing weight variable names
#' @param R A character string containing randomization variable name
#'
#' @return A logical value
#' @keywords internal
#' @examples
#' \dontrun{
#' n <- 200
#' R <- rbinom(n, 1, 0.7)
#' data <- data.frame(
#'   Y = ifelse(R == 1, rnorm(n), NA),
#'   A = ifelse(R == 1, rbinom(n, 1, 0.5), NA),
#'   X1 = rnorm(n),
#'   R = R
#' )
#' check_r_ind(data, Y = "Y", A = "A", X = "X1", W = character(0), R = "R")
#' }
check_r_ind <- function(data,
                        Y,A,X,W,R) {

  # when R=1, Y A X and W should be available. when R=0 something should be missing
  # get rows where R=1
  data_r1 <- data[data[[R]]==1,]

  # get rows where R=0
  data_r0 <- data[data[[R]]==0,]

  # check if Y, A, X and W are all available when R=1
  check_r1 <- all(complete.cases(data_r1[,c(Y,A,X,W)]))

  # check if Y, A, X and W are all missing when R=0
  check_r0 <- all(!complete.cases(data_r0[,c(Y,A,X,W)]))

  # send error messages if r1 and or r0 are not satisfied
  if(!check_r1 & !check_r0) {
    stop("Error: cases where R=1 but Y, A, X and W are not all available and cases where R=0 but Y, A, X and W are not all missing")
  } else if(!check_r1) {
    stop("Error: cases where R=1 but Y, A, X and W are not all available")
  } else if(!check_r0) {
    stop("Error: cases where R=0 but Y, A, X and W are not all missing")
  } else {
    return(TRUE)
  }

}

#' @title Check arguments to drcmd for entry errors
#'
#' @description Checks if the arguments to main function drcmd are correctly specified.
#' Returns TRUE if all checks are passed, throws an error with message otherwise.
#'
#' @param Y A vector or data frame containing outcome values
#' @param A A vector or data frame  containing treatment variable values
#' @param X A data frame containing covariate values
#' @param W A data frame containing proxy variable values
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' n <- 200
#' X <- data.frame(X1 = rnorm(n))
#' W <- data.frame(W1 = rnorm(n))
#' A <- rbinom(n, 1, 0.5)
#' Y <- rnorm(n)
#' check_entry_errors(Y, A, X, W, eem_ind = FALSE, Rprobs = NA, k = 1)
#' }
check_entry_errors <- function(Y,A,X,W,
                               eem_ind,Rprobs,
                               k,cutoff=0.025,cv_folds=5,
                               tml=FALSE,quiet=TRUE,parallel=FALSE,
                               att=FALSE,atc=FALSE) {

 # Make sure Y is a vector
 if (!is.numeric(Y) || !is.null(dim(Y))) {
    stop('Y must be a numeric vector')
  }
  if (any(!is.na(Y) & !is.finite(Y))) {
    stop('Y must contain only finite values or NA')
  }
  if (length(unique(Y[!is.na(Y)])) < 2L) {
    stop('Y must contain at least two distinct observed values')
  }

  # Make sure A is a vector
  if (!is.numeric(A) || !is.null(dim(A))) {
    stop('A must be a numeric vector and 0/1 binary')
  }
  if (any(!is.na(A) & !is.finite(A))) {
    stop('A must contain only finite values or NA')
  }

  # Make sure A is 0/1 binary
  if (!check_binary(A[!is.na(A)])) {
    stop('A must be binary')
  }

  # Make sure X is a data frame
  if (!is.data.frame(X) ) {
    stop('X must be a data frame')
  }
  if (ncol(X) == 0L) {
    stop('X must contain at least one column')
  }

  # Make sure W is a data frame
  if (!is.data.frame(W)  ) {
    stop('W must be a data frame')
  }

  has_nonfinite <- function(data) {
    any(vapply(data, function(x) {
      is.numeric(x) && any(!is.na(x) & !is.finite(x))
    }, logical(1)))
  }
  if (has_nonfinite(X)) {
    stop('Numeric columns in X must contain only finite values or NA')
  }
  if (has_nonfinite(W)) {
    stop('Numeric columns in W must contain only finite values')
  }

  constant_columns <- function(data) {
    names(data)[vapply(data, function(x) {
      length(unique(x[!is.na(x)])) < 2L
    }, logical(1))]
  }
  constant_X <- constant_columns(X)
  constant_W <- constant_columns(W)
  if (length(constant_X) > 0L) {
    stop('Columns in X must contain at least two distinct observed values: ',
         paste(constant_X, collapse = ', '))
  }
  if (length(constant_W) > 0L) {
    stop('Columns in W must contain at least two distinct observed values: ',
         paste(constant_W, collapse = ', '))
  }

  # W contains variables that must be available for every observation
  if (anyNA(W)) {
    stop('W must not contain missing values')
  }

  # Make sure no vars in X are named 'Y' or 'A'
  if (any(colnames(X) %in% c('Y','A','y'))) {
    stop('No variables in X can be named "Y" "y" or "A", which are reserved for outcome and treatment')
  }

  # Make sure no vars in W are named 'Y' or 'A'
  if (any(colnames(W) %in% c('Y','A','y'))) {
    stop('No variables in W can be named "Y" "y" or "A", which are reserved for outcome and treatment')
  }

  # Make sure Y A X and W have same # of observations
  if (length(Y) != nrow(X) | length(Y) != nrow(W) | length(A) != nrow(X) | length(A) != nrow(W) | length(Y) != length(A) | nrow(X) != nrow(W) ) {
    stop('Y, A, X and W must have the same number of observations')
  }

  logical_args <- list(eem_ind=eem_ind,tml=tml,quiet=quiet,
                       parallel=parallel,att=att,atc=atc)
  for (name in names(logical_args)) {
    value <- logical_args[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(name, ' must be TRUE or FALSE')
    }
  }

  Rprobs_missing <- length(Rprobs) == 1L && is.atomic(Rprobs) && is.na(Rprobs)
  if (!Rprobs_missing &&
      (!is.numeric(Rprobs) || length(Rprobs) != length(A) ||
       anyNA(Rprobs) || any(!is.finite(Rprobs)) ||
       any(Rprobs <= 0) || any(Rprobs > 1))) {
    stop('Rprobs must be NA or a vector of probabilities in (0, 1]')
  }

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 1 || k != floor(k) || k > length(Y)) {
    stop('k must be a positive integer no greater than the number of observations')
  }

  max_cv_folds <- if (k == 1L) length(Y) else length(Y) - ceiling(length(Y) / k)
  if (!is.numeric(cv_folds) || length(cv_folds) != 1L ||
      !is.finite(cv_folds) || cv_folds < 2 ||
      cv_folds != floor(cv_folds) || cv_folds > max_cv_folds) {
    stop('cv_folds must be an integer between 2 and the cross-fitting training-sample size')
  }

  if (!is.null(cutoff) &&
      (!is.numeric(cutoff) || length(cutoff) != 1L ||
       !is.finite(cutoff) || cutoff < 0 || cutoff >= 0.5)) {
    stop('cutoff must be NULL or a single number in [0, 0.5)')
  }

  return(TRUE)

}

#' @title Truncate treatment propensity scores
#'
#' @description Truncate propensity scores to interval `[c, 1-c]`
#' @param x A vector of treatment propensity scores
#'
#' @return A vector of treatment propensity scores truncated to interval `[c, 1-c]`
#' @keywords internal
#' @examples
#' \dontrun{
#' x <- c(0.001, 0.3, 0.5, 0.7, 0.999)
#' suppressWarnings(truncate_g(x, cutoff = 0.025))
#' }
truncate_g <- function(x, cutoff=0.025) {
  if (any( (x > 1 - cutoff) | (x < cutoff))) {
    warning(paste0("Propensity scores outside of ", cutoff, " and ", 1-cutoff, ". Truncating to cutoffs"))
  }
  x <- ifelse(x > 1 - cutoff, 1 - cutoff, ifelse(x < cutoff, cutoff, x))
  return(x)
}

#' @title Truncate complete case propensity scores
#'
#' @description Truncate propensity scores to interval `[c, 1-c]`
#' @param x A vector of complete case propensity scores
#'
#' @return A vector of complete propensity scores truncated to interval `[c, 1-c]`
#' @keywords internal
#' @examples
#' \dontrun{
#' x <- c(0.001, 0.3, 0.5, 0.7, 0.999)
#' suppressWarnings(truncate_r(x, cutoff = 0.01))
#' }
truncate_r <- function(x, cutoff=0.01) {
  if (any( (x > 1 - cutoff) | (x < cutoff))) {
    warning(paste0("Complete case probabilities outside of ", cutoff, " and ", 1-cutoff, ". Truncating to cutoffs"))
  }
  x <- ifelse(x > 1 - cutoff, 1 - cutoff, ifelse(x < cutoff, cutoff, x))
  return(x)
}

#' @title Trim vector (for numerical stability)
#'
#' @description Trims values of a vector to avoid numerical instability
#' @param x A vector of values
#' @param val A small value to add to 0 and subtract from 1
#' @return A vector with values trimmed to avoid numerical instability
#' @keywords Internal
trim <- function(x,val=.Machine$double.eps) {
  pmin(pmax(x, val), 1 - val)
}


#' @title Clean SuperLearner libraries
#'
#' @description Internal function for setting SuperLearner libraries before they
#' are passed into estimation procedures. Default libraries from the sl_learners
#' argument are used to fill in missing libraries for any nuisance function with
#' libraries left unspecified
#'
#' @param default_learners Either null, or a character vector containing SuperLearner
#'  libraries to use for estimating all nuisance functions. User can alternatively
#'  specify libraries for each nuisance function for added flexibility
#' @param m_learners Either null, or a character vector containing SuperLearner
#' libraries to be used for the outcome regression
#' @param g_learners Either null, or a character vector containing SuperLearner
#' libraries to be used for the propensity scores
#' @param r_learners Either null, or a character vector containing SuperLearner
#' libraries to be used for the missingness indicator regression
#' @param po_learners Either null, or a character vector containing SuperLearner
#' libraries to be used for the pseudo outcome regression
#'
#' @return A list of
clean_learners <- function(default_learners,
                               m_learners,g_learners,
                               r_learners,po_learners) {

  out_list <- list(default_learners=default_learners,
                   m_learners=m_learners,
                   g_learners=g_learners,
                   r_learners=r_learners,
                   po_learners=po_learners)
  if (!is.null(default_learners)) {
    out_list <- lapply(out_list, function(x) if (is.null(x)) default_learners else x)
  } else {
    if( any(sapply(list(m_learners,g_learners,r_learners,po_learners),
                   is.null)) ) {
      stop('Missing learner specifications. Make sure all learners are specified, either through specifying default values through default_learners or specifying learners for each nuisance function')
    }
  }
  return(out_list)
}

#' @title check_binary
#' @description Check if the outcome variable is 0/1 binary. Return 1 if true, 0
#' if false
#' @param x A numeric vector
#' @return A logical value
check_binary <- function(x) {
  x <- x[!is.na(x)]
  length(x) > 0 && all(x %in% c(0, 1))
}

# Fractional observation weights do not represent binomial trial counts.
# An all-zero ensemble, however, is not an acceptable probability model.
with_probability_fit_checks <- function(expr, model) {
  withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (msg %in% c("All algorithms have zero weight",
                   "All metalearner coefficients are zero, predictions will all be equal to 0")) {
      stop(model, " produced an all-zero SuperLearner ensemble; ",
           "revise the learner library or the cross-validation specification. ",
           "Probability truncation cannot repair this fit.", call. = FALSE)
    }
    if (identical(msg, "non-integer #successes in a binomial glm!")) {
      invokeRestart("muffleWarning")
    }
  })
}

# Consolidate SuperLearner's repeated zero-coefficient warnings into one
# diagnostic per augmentation fit. Other fitting warnings remain visible.
fit_augmentation_sl <- function(...) {
  fit <- withCallingHandlers(SuperLearner::SuperLearner(...), warning = function(w) {
    if (conditionMessage(w) %in% c("All algorithms have zero weight",
        "All metalearner coefficients are zero, predictions will all be equal to 0")) {
      invokeRestart("muffleWarning")
    }
  })
  if (!length(fit$coef) || any(!is.finite(fit$coef))) {
    stop("Augmentation regression produced invalid ensemble coefficients", call. = FALSE)
  }
  if (all(fit$coef == 0)) {
    warning(structure(list(
      message = paste("Augmentation regression selected an all-zero ensemble;",
                      "its predictions are zero. Consider revising the learner library."),
      call = NULL), class = c("drcmd_zero_augmentation", "warning", "condition")))
  }
  fit
}

predict_augmentation_sl <- function(fit, newdata) {
  if (all(fit$coef == 0)) return(rep(0, nrow(newdata)))
  pred <- predict(fit, newdata = newdata)$pred
  if (any(!is.finite(pred))) {
    stop("Augmentation regression produced nonfinite predictions", call. = FALSE)
  }
  pred
}

check_probability_predictions <- function(x, model) {
  if (!length(x) || any(!is.finite(x) | x < 0 | x > 1)) {
    stop(model, " produced invalid probability predictions", call. = FALSE)
  }
  if (all(x == 0) || all(x == 1)) {
    stop(model, " produced only boundary probabilities; ",
         "revise the learner library or the cross-validation specification. ",
         "Probability truncation cannot repair this fit.", call. = FALSE)
  }
  invisible(x)
}

binary_cv_control <- function(y, weights, V, model) {
  shuffle <- function(x) if (length(x) > 1L) sample(x) else x
  rows <- which(is.finite(weights) & weights > 0)
  counts <- tabulate(y[rows] + 1L, nbins = 2L)
  if (any(counts < 2L)) {
    stop(model, " requires at least two positive-weight training observations ",
         "in each class for cross-validation", call. = FALSE)
  }

  valid_rows <- vector("list", V)
  for (value in 0:1) {
    class_rows <- shuffle(rows[y[rows] == value])
    fold <- rep(seq_len(V), length.out = length(class_rows))
    for (i in seq_len(V)) {
      valid_rows[[i]] <- c(valid_rows[[i]], class_rows[fold == i])
    }
  }

  remaining <- shuffle(setdiff(seq_along(y), rows))
  fold <- rep(seq_len(V), length.out = length(remaining))
  for (i in seq_len(V)) {
    valid_rows[[i]] <- shuffle(c(valid_rows[[i]], remaining[fold == i]))
  }

  list(V = V, validRows = valid_rows)
}


#' @title Create folds for cross-fitting
#' @description Given length of data, creates k folds and returns a list of all
#' possible train-test pairs
#'
#' @param n Length of data
#' @param k Number of desired folds
#'
#' @return A list of train-test pairs
create_folds <- function(n, k) {

  if (k==1) {
    return(list(test=1:n, train=1:n))
  }

  # Randomly split the data into K pieces
  indices <- sample(seq_len(n))
  split_indices <- split(indices, cut(seq_along(indices), breaks = k, labels = FALSE))

  # Create train-test pairs
  splits <- lapply(seq_len(k), function(i) {
    test_indices <- split_indices[[i]]
    train_indices <- setdiff(indices, test_indices)
    list(test = test_indices, train = train_indices)
  })

  return(splits)
}

#' @title List SuperLearner libraries
#'
#' @description List all available SuperLearner libraries
#' @return A character vector of all available SuperLearner libraries
#' @export
#' @examples
#' \dontrun{
#' libs <- get_sl_libraries()
#' head(libs)
#' }
get_sl_libraries <- function() {
  sink(tempfile())  # Redirect console output
  all_wrappers <- suppressMessages(SuperLearner::listWrappers())
  sink()  # Restore console output
  SL_wrappers <- c(all_wrappers[grep("^SL\\.", all_wrappers)], "SL.hal9001")

  ''
  return(SL_wrappers)
}


get_clean_context <- function(calls) {
  for (i in rev(seq_along(calls))) {
    fn <- calls[[i]]
    if (is.call(fn) && is.symbol(fn[[1]])) {
      fname <- as.character(fn[[1]])
      if (!any(grepl(fname, c("withCallingHandlers", "eval", "doWithOneRestart",
                              "withRestarts", "signalCondition", "structure",
                              "base::quote", "get_clean_context", "muffleWarning")))) {
        return(deparse1(fn))
      }
    }
  }
  return("unknown")
}
