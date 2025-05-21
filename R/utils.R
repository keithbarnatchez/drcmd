# helper functions for the drcmd function

#' @title find_missing_pattern
#' @description Find the missing pattern in the data
#' @param data A data frame
#' @param Y A vector or data frame containing outcome values
#' @param A A vector or data frame  containing treatment variable values
#' @param X A data frame containing covariate values
#' @param W A data frame containing proxy variable values
#' @param R A vector containing missingness indicator variable
#'
#' @return A character string containing the missing pattern
#' @export
find_missing_pattern <- function(Y,A,X,W) {

  # Combine variables into a single data frame
  data <- cbind(X,W,
                data.frame(Y=Y,A=A))

  # find variables that are never missing
  never_missing <- colnames(data)[apply(data,2,function(x) sum(is.na(x))==0)]
  if(length(never_missing) == 0) {
    stop("Error: drcmd requires data to have at least one variable that is never missing")
  }

  # find variables that are sometimes missing
  sometimes_missing <- colnames(data)[apply(data,2,function(x) sum(is.na(x))>0)]

  # make variable that indicates complete cases and df of all complete data
  R <- as.numeric(apply(data,1,function(x) sum(is.na(x))==0))
  if (all(R == 0)) {
    stop("Error: drcmd requires data to have at least one complete case")
  }
  if (mean(R) < 0.01) {
    warning('Small number of complete cases. Results may be unstable')
  }

  Z <- data[,colnames(data)[apply(data,2,function(x) sum(is.na(x))==0)],drop=FALSE]

  # If 'Y' is in Z rename to 'y' (superlearner doesn't allow covariates named Y)
  if('Y' %in% colnames(Z)) {
    colnames(Z)[colnames(Z) == 'Y'] <- 'y'
  }

  # If any values of Y, X or Z equal NA, set them to 0
  # Doing this since some learners don't support NAs (these obs get 0 weight anyway)
  X[is.na(X)] <- 0
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
#' @export
check_r_ind <- function(data,
                        Y,A,X,W,R) {

  # when R=1, Y A X and W should be available. when R=0 something should be missing
  # get rows where R=1
  data_r1 <- data[data[[R]]==1,]

  # get rows where R=0
  data_r0 <- data[data[[R]]==0,]

  # check if Y, A, X and W are all available when R=1
  check_r1 <- complete.cases(data_r1[,c(Y,A,X,W)])

  # check if Y, A, X and W are all missing when R=0
  check_r0 <- !complete.cases(data_r0[,c(Y,A,X,W)])

  # send error messages if r1 and or r0 are not satisfied
  if(!check_r1) {
    stop("Error: cases where R=1 but Y, A, X and W are not all available")
  } else if(!check_r0) {
    stop("Error: cases where R=0 but Y, A, X and W are not all missing")
  } else if(!check_r1 & !check_r0) {
    stop("Error: cases where R=1 but Y, A, X and W are not all available and cases where R=0 but Y, A, X and W are not all missing")
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
#' @param R A vector containing missingness indicator variable
#'
#' @export
check_entry_errors <- function(Y,A,X,W,R,
                               eem_ind,Rprobs,
                               k,nboot) {

 # Make sure Y is a vector
  if (!is.double(Y) & !is.integer(Y)) {
    stop('Y must be a numeric vector')
  }

  # Make sure A is a vector
  if (!is.vector(A) & !is.integer(A) ) {
    stop('A must be a numeric vector and 0/1 binary')
  }

  # Make sure A is 0/1 binary
  if (!check_binary(A[!is.na(A)])) {
    stop('A must be binary')
  }

  # Make sure X is a data frame
  if (!is.data.frame(X) ) {
    stop('X must be a data frame')
  }

  # Make sure W is a data frame
  if (!is.data.frame(W)  ) {
    stop('W must be a data frame')
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

  # check that eem_ind is a logical
  if (!is.logical(eem_ind) | (length(eem_ind)>1) ) {
    stop('eem_ind must be a logical')
  }

  # check that Rprobs is a vector of values between 0 and 1 inclusive
  if (!any(is.na(Rprobs))) {
    if (!is.numeric(Rprobs) | any(Rprobs < 0) | any(Rprobs > 1)  | length(Rprobs)!=length(A) ) {
      stop('When specicified, Rprobs must be a vector of sampling probabilities between 0 and 1 inclusive')
    }
  }
  # make sure cross-fitting folds is an integer
  if (!is.numeric(k)) {
    stop('k must be a an integer')
  } else{
    if (floor(k) != k) {
      stop('k must be an integer')
    }
  }

  if (k > length(Y)) {
    stop('k must be less than the number of observations')
  }

  if (!is.numeric(nboot) | nboot<0) {
    stop('nboot must be a an integer')
  } else{
    if (floor(nboot) != nboot) {
      stop('nboot must be an integer')
    }
  }

  return(TRUE)

}

#' @title Truncate treatment propensity scores
#'
#' @description Truncate propensity scores to interval [c, 1-c]
#' @param x A vector of treatment propensity scores
#'
#' @return A vector of treatment propensity scores truncated to interval [c, 1-c]
#' @export
#'
truncate_g <- function(x, cutoff=0.025) {
  if (any( (x > 1 - cutoff) | (x < cutoff))) {
    warning(cat("Propensity scores outside of ",cutoff," and ",1-cutoff,". Truncating to cutoffs\n"))
  }
  x <- ifelse(x > 1 - cutoff, 1 - cutoff, ifelse(x < cutoff, cutoff, x))
  return(x)
}

#' @title Truncate complete case propensity scores
#'
#' @description Truncate propensity scores to interval [c, 1-c]
#' @param x A vector of complete case propensity scores
#'
#' @return A vector of complete propensity scores truncated to interval [c, 1-c]
#' @export
#'
truncate_r <- function(x, cutoff=0.01) {
  if (any( (x > 1 - cutoff) | (x < cutoff))) {
    warning(cat("Complete case probabilities outside of ",cutoff," and ",1-cutoff,". Truncating to cutoffs\n"))
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
trim <- function(x,val=.Machine$double.neg.eps) {
  x[x==0] <- x+val
  x[x==1] <- x-val
  return(x)
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
#' @param m_sl_learners Either null, or a character vector containing SuperLearner
#' libraries to be used for the outcome regression
#' @param g_sl_learners Either null, or a character vector containing SuperLearner
#' libraries to be used for the propensity scores
#' @param r_sl_learners Either null, or a character vector containing SuperLearner
#' libraries to be used for the missingness indicator regression
#' @param po_sl_learners Either null, or a character vector containing SuperLearner
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
  # check if x is binary
  if(all(x %in% c(0,1))) {
    return(TRUE)
  } else {
    if (min(x) == 0 & max(x) == 1) { # if y normed to unit interval, treat as binary
      return(TRUE)
    }
    return(FALSE)
  }
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
get_sl_libraries <- function() {
  sink(tempfile())  # Redirect console output
  all_wrappers <- suppressMessages(SuperLearner::listWrappers())
  sink()  # Restore console output
  SL_wrappers <- c(all_wrappers[grep("^SL\\.", all_wrappers)], "SL.hal9001")

  ''
  return(SL_wrappers)
}

#' @title Clean nuisance function output from crossfit procedure
#'
#' @description Transforms nuisance output across folds into a dataframe of avgerge
#' prediction at each point for each nuisance function
#'
#' @param results Nuisance functions output from drcmd_est
#' @return A dataframe of averaged nuisance estimates across all folds
#' @export
clean_crossfit_nuis <- function(results) {
  1
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

