# helper functions for the drcmd function

#' @title find_missing_pattern
#' @description Find the missing pattern in the data
#' @param data A data frame
#' @param Y A character string containing outcome variable name
#' @param A A character string containing treatment variable name
#' @param X A character vector containing covariate variable names
#' @param W A character vector containing weight variable names
#' @param R A character string containing randomization variable name
#'
#' @return A character string containing the missing pattern
#' @export
find_missing_pattern <- function(data,
                                 Y,A,X,W,R) {

  # find variables that are never missing
  never_missing <- colnames(data)[apply(data,2,function(x) sum(is.na(x))==0)]
  Rstr <- R

  if (is.na(R)) {
    # make variable that indicates complete cases
    Rvals <- as.numeric(apply(data,1,function(x) sum(is.na(x))==0))
    never_missing <- colnames(data)[apply(data,2,function(x) sum(is.na(x))==0)]
    Rstr <- 'R'
  } else {
    Rvals <- data[,R]
    # remove R from never_missing
    never_missing <- never_missing[!never_missing %in% R]
  }

  return(list(Z=never_missing,
              Rstr=Rstr,
              Rvals=Rvals))

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

check_entry_errors <- function(Y,A,X,W,R, hal_ind,sl.lib,eem_ind,Rprobs,k) {

  # Make sure Y,A and R are characters
  if ( !is.character(Y) | (length(Y)>1) ) {
    stop('Y must be a character string')
  }
  if (!is.character(A) | (length(A)>1) ) {
    stop('A must be a string')
  }
  if (is.character(R) & (length(R)>1) ) {
    stop('R must be a string or not provided')
  }

  # Make sure X and W are character vectors
  if (!is.character(X) ) {
    stop('X must be a character vector')
  }
  if (!is.character(W) ) {
    stop('W must be a character vector')
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


#'
#'
#'
#'
#'
#'
#'
#'
create_folds <- function(data, k) {

  n <- nrow(data)
  if (k==1) {
    return(list(test=1:n, train=1:n))
  }

  # Step 1: Randomly split the data into K pieces
  set.seed(123)  # Optional: set a seed for reproducibility
  indices <- sample(seq_len(n))
  split_indices <- split(indices, cut(seq_along(indices), breaks = k, labels = FALSE))

  # Step 2: Create train-test pairs
  splits <- lapply(seq_len(k), function(i) {
    test_indices <- split_indices[[i]]
    train_indices <- setdiff(indices, test_indices)
    list(test = test_indices, train = train_indices)
  })

  return(splits)
}


