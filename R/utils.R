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
  data <- cbind(X,W,data.frame(Y=Y,A=A))

  # find variables that are never missing
  never_missing <- colnames(data)[apply(data,2,function(x) sum(is.na(x))==0)]

  # find variables that are sometimes missing
  sometimes_missing <- colnames(data)[apply(data,2,function(x) sum(is.na(x))>0)]

  # make variable that indicates complete cases and df of all complete data
  R <- as.numeric(apply(data,1,function(x) sum(is.na(x))==0))
  Z <- data[,colnames(data)[apply(data,2,function(x) sum(is.na(x))==0)]]

  # If any values of Y, X or Z equal NA, set them to 0
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

check_entry_errors <- function(Y,A,X,W,R,
                               hal_ind,sl.lib,eem_ind,Rprobs,k) {

 # Make sure Y is a vector
  if (!is.vector(Y) | (length(Y)>1) ) {
    stop('Y must be a vector')
  }

  # Make sure A is a vector
  if (!is.vector(A) | (length(A)>1) ) {
    stop('A must be a vector')
  }

  # Make sure A is 0/1 binary
  if (!check_binary(A)) {
    stop('A must be binary')
  }

  # Make sure X is a data frame
  if (!is.data.frame(X) | (length(X)>1) ) {
    stop('X must be a data frame')
  }

  # Make sure W is a data frame
  if (!is.data.frame(W) | (length(W)>1) ) {
    stop('W must be a data frame')
  }

  # Check that hal_ind is a logical
  if (!is.logical(hal_ind) | (length(hal_ind)>1) ) {
    stop('hal_ind must be a logical')
  } else if ( (hal_ind == FALSE) & is.na(sl.lib) ) {
    stop('If not using HAL, must provide superlearner library')
  }

  # check that eem_ind is a logical
  if (!is.logical(eem_ind) | (length(eem_ind)>1) ) {
    stop('eem_ind must be a logical')
  }

  return(TRUE)

}


#' @title Clean SuperLearner libraries
#'
#' @description Internal function for setting SuperLearner libraries before they
#' are passed into estimation procedures. Default libraries from the sl_learners
#' argument are used to fill in missing libraries for any nuisance function with
#' libraries left unspecified
#'
#' @param sl_learners Either null, or a character vector containing SuperLearner
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
clean_sl_libraries <- function(sl_learners,
                               m_sl_learners,g_sl_learners,
                               r_sl_learners,po_sl_learners) {

  out_list <- list(sl_learners=sl_learners,
                   m_sl_learners=sl_learners,
                   g_sl_learners=sl_learners,
                   r_sl_learners=sl_learners,
                   po_sl_learners=sl_learners)
  if (!is.null(sl_learners)) {
    out_list <- lapply(out_list, function(x) if (is.null(x)) sl_learners else x)
  } else {
    if( any(sapply(list(m_sl_learners,g_sl_learners,r_sl_learners,po_sl_learners),
                   is.null)) ) {
      stop('Missing superlearner libraries. Make sure all libraries are specified, either through specifying default values through sl_learners or specifying learners for each nuisance function')
    }
  }

  return(out_list)

}

#' @title check_binary
#' @description Check if the outcome variable is 0/1 binary. Return 1 if true, 0
#' if false
#' @param x A numeric vector
#' @return A logical value
#' @export
check_binary <- function(x) {
  # check if x is binary
  if(all(x %in% c(0,1))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

make_output_obj <- function() {

}


#' @title Create folds for cross-fitting
#' @description Given length of data, creates k folds and returns a list of all
#' possible train-test pairs
#'
#' @param n Length of data
#' @param k Number of desired folds
#'
#' @return A list of train-test pairs
#' @export
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



#' @title Clean output from crossfit procedure
#'
#' @description Transforms crossfit output into estimates with SEs for each
#' estimand
#'
#' @param results Output from crossfit procedure
#' @return A list of estimates with SEs
#' @export
clean_crossfit_output <- function(results) {

}


