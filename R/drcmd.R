#' @title drcmd
#' @description Main function for the drcmd package.
#' @param data A data frame
#' @param Y A character string containing outcome variable name
#' @param A A character string containing treatment variable name
#' @param X A character vector containing covariate names
#' @param W A character vector containing variable names solely used for imputation
#' @param R (optional) A character string specifying the missingness indicator, where 0 indicates missing data. If not specified, the function will identify the missingness indicator
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl.lib A character string indicating the superlearner library to use
#' @param eem_ind A logical indicating whether to use ensemble of estimators
#' @param Rprobs A vector of probabilities for R
#' @param k A numeric indicating the number of folds for cross-fitting
#' @param nboot A numeric indicating the number of desired bootstrap samples. If >0, uses bootstrap to obtain SEs. If =0, uses asymptotic analytical SEs.
#' @return A list of results
#'
#' @export
drcmd <- function(data, Y=NA, A=NA, X=NA, W=NA, R=NA,
                  hal_ind=TRUE, sl.lib=NA, eem_ind=FALSE, Rprobs=NA, k=1,
                  nboot=0) {

  require(SuperLearner)

  # Throw errors if anything is entered incorrectly
  check_entry_errors(Y,A,X,W,R, hal_ind,sl.lib,eem_ind,Rprobs)

  # Identify missing data structure
  V <- find_missing_pattern(data,Y,A,X,W,R)
  Z <- V$Z
  R <- V$Rstr
  data[,R] <- V$Rvals

  res <- drcmd_est(data,Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs,k)

  # Return results
  return(res)

}

#' @title drcmd_est
#'
#' @description Outer function for obtaining point estimates and standard errors
#' @param data A data frame
#' @param Y A character string
#' @param A A character string
#' @param X A character vector
#' @param Z A character vector
#' @param R A character vector
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl.lib A character string indicating the superlearner library to use
#' @param eem_ind A logical indicating whether to use ensemble of estimators
#' @param Rprobs A vector of probabilities for R
#' @param k A numeric indicating the number of folds for cross-fitting
#'
#' @return A list of results
#' @export
drcmd_est <- function(data,Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs,k) {

  # Divide the data into k random "train-test" splits
  if (k>1) {
    splits <- create_folds(data,k)

    # ests <- drcmd_est_fold(data,splits[[1]],Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs)

    # Use lapply to get point ests for each fold
    ests <- lapply(1:k, function(i) drcmd_est_fold(data,splits[[i]],Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs))
    ests <- colMeans(do.call(rbind, lapply(ests, as.data.frame)))
    return(ests)
  } else {
    splits <- create_folds(data,k)
    ests <- drcmd_est_fold(data,splits,Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs)
    return(ests)
  }
}

#' @title drcmd_est_fold
#'
#' @description Outer function for obtaining point estimates for single fold
#' @param data A data frame
#' @param splits A list of train/test splits
#' @param Y A character string
#' @param A A character string
#' @param X A character vector
#' @param Z A character vector
#' @param R A character vector
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl.lib A character string indicating the superlearner library to use
#' @param eem_ind A logical indicating whether to use ensemble of estimators
#' @param Rprobs A vector of probabilities for R
#' @param k A numeric indicating the number of folds for cross-fitting
#'
#' @return A list of results
#' @export
drcmd_est_fold <- function(data,splits,Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs) {



  # Get the training and test data
  train <- data[splits$train,]
  test <- data[splits$test,]

  # Get the nuisance estimates
  nuisance_ests <- get_nuisance_ests(train,Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs)

  # Form phi
  phi_hat <- get_phi_hat(test,Y,A,X,Z,R,nuisance_ests$g_hat,nuisance_ests$m_a_hat,
                      nuisance_ests$kappa_hat,hal_ind,sl.lib)

  # Get varphi
  phi_1_hat <- phi_hat$phi_1_hat ; phi_0_hat <- phi_hat$phi_0_hat
  varphi_hat <- est_varphi(test,R,Z,phi_1_hat,phi_0_hat,hal_ind,sl.lib)

  # Form final est for this fold
  ests <- est_psi(test, R, Z, nuisance_ests$kappa_hat, phi_hat,varphi_hat)

  # Get the standard error

  # Return results
  return(ests)

}


#' @title get_nuisance_ests
#'
#' @description Function for obtaining estimates of all relevant nuisance functions
#' @param data A data frame
#' @param Y A character string containing outcome variable name
#' @param A A character string containing treatment variable name
#' @param X A character vector containing covariate names
#' @param Z A character vector containing first stage variable names
#' @param R A character vector containing missing data indicator
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl.lib A character string indicating the superlearner library to use
#' @param eem_ind A logical indicating whether to use ensemble of estimators
#' @param Rprobs A vector of probabilities for R
#'
#' @return A list of nuisance estimates
#' @export
get_nuisance_ests <- function(data,Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs) {

  kappa_hat <- est_kappa(data,Z,R,hal_ind,sl.lib)

  m_a_hat <- est_m_a(data,Y,A,X,R,kappa_hat,hal_ind,sl.lib)
  g_hat <- est_g(data,A,X,R,kappa_hat,hal_ind,sl.lib)

  return(list(kappa_hat=kappa_hat,m_a_hat=m_a_hat,g_hat=g_hat))

}

#' @title get_phi_hat
#'
#' @description Function for obtaining estimate of E[phi|Z], where phi is the
#' full data EIF
#'
#' @param data A data frame
#' @param Y A character string containing outcome variable name
#' @param A A character string containing treatment variable name
#' @param X A character vector containing covariate variable names
#' @param Z A character vector containing the names of the variables in Z
#' @param R A character string containing randomization variable name
#' @param g_hat A fitted propensity score model
#' @param m_a_hat A fitted outcome model
#' @param kappa_hat A fitted sampling probs model
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl.lib A character string indicating the superlearner library to use
#' @return A list containing the estimated EIFs for A=1 and A=0
#'
#' @export
#'
get_phi_hat <- function(data, Y, A, X, Z, R, g_hat, m_a_hat, kappa_hat,
                        hal_ind,sl.lib) {


  orig_A <- data[,A]
  # First, get predicted values for current data
  if (hal_ind==TRUE) { # estimate via HAL
    # get fitted outcome values under A=1 and A=0
    data[,A] <- 1 ; m_1_hat <- predict(m_a_hat, newdata=data)
    data[,A] <- 0 ; m_0_hat <- predict(m_a_hat, newdata=data)

    # get propensity scores
    g_hat <- predict(g_hat, newdata=data)
  } else { # estimate via SL *** need dynamic family ***
    data[,A] <- 1 ; m_1_hat <- predict(m_a_hat, newdata=data)$pred
    data[,A] <- 0 ; m_0_hat <- predict(m_a_hat, newdata=data)$pred

    # get propensity scores
    g_hat <- predict(g_hat, newdata=data)$pred
  }
  data[,A] <- orig_A

  # get fitted values under A=1
  phi_1_hat <- m_1_hat + data[,A]*(data[,Y] - m_1_hat)/g_hat
  phi_0_hat <- m_0_hat + (1-data[,A])*(data[,Y] - m_0_hat)/(1-g_hat)

  # Set values where R==0 to 0
  phi_1_hat[data[,R]==0] <- 0
  phi_0_hat[data[,R]==0] <- 0

  return(list(phi_1_hat=phi_1_hat,phi_0_hat=phi_0_hat))
}

#' @title est_varphi
#' @description Function for obtaining estimate of E[phi|Z]
#' @param data A data frame
#' @param R A character string containing randomization variable name
#' @param phi_1_hat A numeric vector containing the fitted values of phi under A=1
#' @param phi_0_hat A numeric vector containing the fitted values of phi under A=0
#' @param hal_ind A logical value indicating whether to estimate via HAL
#' @param sl.lib A character vector containing the names of the superlearner algorithms
#' @return A list containing the estimate of E[phi|Z]
#' @export
est_varphi <- function(data, R, Z,
                       phi_1_hat, phi_0_hat,
                       hal_ind,
                       sl.lib) {

  if (hal_ind==TRUE) { # estimate via HAL
    varphi_1_hat <- fit_hal(Y=phi_1_hat,X=data[,Z],weights=data[,R])
    varphi_0_hat <- hal9001::fit_hal(Y=phi_0_hat,X=data[,Z])
  } else { # estimate via SL
    varphi_1_hat <- SuperLearner::SuperLearner(Y=phi_1_hat,X=data[,Z],
                                               family=gaussian(),SL.library=sl.lib,
                                               obsWeights=data[,R])
    varphi_0_hat <- SuperLearner::SuperLearner(Y=phi_0_hat,X=data[,Z],
                                               family=gaussian(),SL.library=sl.lib,
                                               obsWeights=data[,R])
  }

  return(list(varphi_1_hat=varphi_1_hat,varphi_0_hat=varphi_0_hat))
}

#' @title est_psi
#' @description Function for obtaining estimate of the full data EIC, psi
#' @param data A data frame
#' @param R A character string containing randomization variable name
#' @param Z A character vector containing the names of the variables in Z
#' @param kappa_hat A numeric vector containing the fitted values of kappa
#' @param phi_hat A list containing the estimate of E[phi|Z]
#' @param varphi_hat A list containing the estimate of E[phi|Z]
#'
#' @return A numeric value containing the estimate of E[psi]
#' @export
est_psi <- function(data, R, Z,
                    kappa_hat, phi_hat,varphi_hat) {

  if (hal_ind) {
    kappa_hat <- predict(kappa_hat, newdata=data)
  } else {
    kappa_hat <- predict(kappa_hat, newdata=data)$pred
  }

  phi_1_hat <- phi_hat$phi_1_hat
  phi_0_hat <- phi_hat$phi_0_hat
  varphi_1_hat <- predict(varphi_hat$varphi_1_hat, newdata=data)$pred
  varphi_0_hat <- predict(varphi_hat$varphi_0_hat, newdata=data)$pred

  psi_1_hat <- mean((data[,R]/kappa_hat)*phi_1_hat - (data[,R]/kappa_hat - 1)*varphi_1_hat)
  psi_0_hat <- mean((data[,R]/kappa_hat)*phi_0_hat - (data[,R]/kappa_hat - 1)*varphi_0_hat)

  psi_hat_ate <- psi_1_hat - psi_0_hat

  return(list(psi_1_hat=psi_1_hat,psi_0_hat=psi_0_hat,psi_hat_ate=psi_hat_ate))

}



#' @title est_m_a
#'
#' @description Function for obtaining estimate of E[Y|A=a,X]
#' @param data A data frame
#' @param Y A character string containing outcome variable name
#' @param A A character string containing treatment variable name
#' @param X A character vector containing the names of the variables in X
#' @param R A character string containing randomization variable name
#' @param kappa_hat A numeric vector containing the fitted values of kappa
#' @param hal_ind A logical value indicating whether to estimate via HAL
#' @param sl.lib A character vector containing the names of the superlearner algorithms
#' @return A list containing the estimate of E[Y|A=a,X]
#' @export
est_m_a <- function(data, Y, A, X, R,
                    kappa_hat,
                    hal_ind,
                    sl.lib) {

  # if Y ever equals NA in data, set those vals to 0
  data[is.na(data[,Y]),Y] <- 0

  # note: need to
  yfam <- check_binary(data[,Y])

  if (hal_ind==TRUE) { # estimate via HAL
    kappa_hat <- predict(kappa_hat, newdata=data)
    m_a_hat<- hal9001::fit_hal(Y=data[,Y],X=data[,c(X,A)],
                               weights=data[,R]/kappa_hat)

     # get fitted values under A=1
    data[,A] <- 1
    m_1_hat <- predict(m_a_hat, newdata=data)
    data[,A] <- 0
    m_0_hat <- predict(m_a_hat, newdata=data)
  } else { # estimate via SL *** need dynamic family ***
    kappa_hat <- predict(kappa_hat, newdata=data)$pred
    m_a_hat <- SuperLearner::SuperLearner(Y=data[,Y],X=data[,c(X,A)],
                                          SL.library=sl.lib,
                                          obsWeights=data[,R]/kappa_hat)
    # get fitted values under A=1
    data[,A] <- 1
    m_1_hat <- predict(m_a_hat, newdata=data)$pred
    data[,A] <- 0
    m_0_hat <- predict(m_a_hat, newdata=data)$pred
  }

  return(m_a_hat)
  # return(list(m_1_hat=m_1_hat,m_0_hat=m_0_hat,m_a_hat=m_a_hat))

}


#' @title est_g
#'
#' @description Function for obtaining estimate of E[Y|A=a,X,W]
#' @param data A data frame
#' @param Y A character string containing outcome variable name
#' @param A A character string containing treatment variable name
#' @param X A character vector containing covariate variable names
#' @param W A character vector containing weight variable names
#' @param R A character string containing randomization variable name
#' @param hal_ind A logical value indicating whether to estimate via HAL
#' @param sl.lib A character vector containing the names of the superlearner algorithms
#'
#' @return A list containing the estimate of E[Y|A=a,X,W]
#' @export
est_g <- function(data, A, X, R, kappa_hat,
                   hal_ind,
                   sl.lib) {

  if (hal_ind==TRUE) { # estimate via HAL
    kappa_hat <- predict(kappa_hat, newdata=data)
    g_hat<- hal9001::fit_hal(Y=data[,A],X=data[,X],weights=data[,R]/kappa_hat)
    # get fitted values
    g_hats <- predict(g_hat)
  } else { # estimate via SL
    kappa_hat <- predict(kappa_hat, newdata=data)$pred
    g_hat <- SuperLearner::SuperLearner(Y=data[,A],X=data.frame(data[,X,drop=FALSE]),
                                        family=binomial(),SL.library=sl.lib,
                                        obsWeights=data[,R]/kappa_hat)
    # get fitted values
    g_hats <- predict(g_hat)$pred
  }
  return(g_hat)
  # return(list(g_hat=g_hat,g_hats=g_hats))
}

# note: handle user supplied r probs outside of this function
est_kappa <- function (data, Z, R,
                       hal_ind,
                       sl.lib) {

  # if estimating via highly adaptive lasso
  if (hal_ind==TRUE) { # estimate via HAL
    kappa_hat <- hal9001::fit_hal(Y=data[,R],X=data[,Z], family = 'binomial')
    # kappa_hat <- predict(kappa_hat, new_data=data)
  }

  # if estimating via superlearner
  if (hal_ind==FALSE) { # estimate via SL
    loadNamespace("SuperLearner")
    kappa_hat <- SuperLearner::SuperLearner(Y=data[,R],X=data[,Z],
                                            family=binomial(),SL.library='SL.glm')
    # kappa_hat <- predict(kappa_hat, newdata=data)$pred

  }

  return(kappa_hat)

}


#' @title est_varphi
#'
#' @description Function for obtaining estimate of E[phi|Z]
#' @param data A data frame
#' @param Y A character string containing outcome variable name
#' @param A A character string containing treatment variable name
#' @param X A character vector containing covariate variable names
#' @param Z A character vector containing the names of the variables in Z
#' @param R A character string containing randomization variable name
#' @param phi_hat A list containing the estimate of phi
#' @param eem_ind A logical value indicating whether to estimate via EEM
#'
#' @return A list containing the estimate of E[phi|Z]
#' @export
est_varphi_main <- function(data, R,
                       phi_1_hat, phi_0_hat,
                       eem_ind,
                       hal_ind,
                       sl.lib){



  if (eem_ind==TRUE) { # estimate via EEM
    return(est_varphi_eem(data, R, Z,
                          phi_1_hat, phi_0_hat,
                          hal_ind,
                          sl.lib))
  } else {
    return(est_varphi(data, R, Z,
                          phi_1_hat, phi_0_hat,
                          hal_ind,
                          sl.lib))
  }


}

