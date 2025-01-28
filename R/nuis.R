#' @title Obtain nuisance function estimates
#'
#' @description Function for obtaining estimates of all relevant nuisance functions
#' @param idx Indices to carry out estimation over
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param Z Dataframe containing all non-missing variables
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl_learners A character string indicating the superlearner library to use
#' @param eem_ind A logical indicating whether to use ensemble of estimators
#' @param Rprobs A vector of probabilities for R
#'
#' @return A list of nuisance estimates
#' @export
get_nuisance_ests <- function(idx,Y,A,X,Z,R,hal_ind,sl_learners,Rprobs) {

  kappa_hat <- Rprobs
  if (is.na(Rprobs)) {
    kappa_hat <- est_kappa(idx,Z,R,hal_ind,sl_learners)
  }

  m_a_hat <- est_m_a(idx,Y,A,X,R,kappa_hat,hal_ind,sl_learners)
  g_hat <- est_g(idx,A,X,R,kappa_hat,hal_ind,sl_learners)

  return(list(kappa_hat=kappa_hat,m_a_hat=m_a_hat,g_hat=g_hat))

}


#' @title Function for estimating outcome regression
#'
#' @description Function for obtaining estimate of E[Y|A=a,X]
#' @param idx Indices to carry out estimation over
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param kappa_hat A numeric vector containing the fitted values of kappa
#' @param hal_ind A logical value indicating whether to estimate via HAL
#' @param sl_learners A character vector containing the names of the superlearner algorithms
#' @return A list containing the estimate of E[Y|A=a,X]
#' @export
est_m_a <- function(idx, Y, A, X, R,
                    kappa_hat,
                    hal_ind,
                    sl_learners) {

  # note: need to
  yfam <- check_binary(Y)

  data <- cbind(X,A)
  if (hal_ind==TRUE) { # estimate via HAL
    m_a_hat<- hal9001::fit_hal(Y=Y[idx],X=data[idx,],
                               weights=R[idx]/kappa_hat[idx])
    # get fitted values under A=1
    data[,ncol(data)] <- 1
    m_1_hat <- predict(m_a_hat, new_data=data)
    data[,ncol(data)] <- 0
    m_0_hat <- predict(m_a_hat, new_data=data)
  } else { # estimate via SL *** need dynamic family ***
    m_a_hat <- SuperLearner::SuperLearner(Y=Y[idx],X=data[idx,],
                                          SL.library=sl_learners,
                                          obsWeights=R[idx]/kappa_hat[idx])
    # get fitted values under A=1
    data[,ncol(data)] <- 1
    m_1_hat <- predict(m_a_hat, newdata=data)$pred
    data[,ncol(data)] <- 0
    m_0_hat <- predict(m_a_hat, newdata=data)$pred
  }

  return(m_a_hat)
}


#' @title Estimate the propensity score
#'
#' @description Function for obtaining propensity score estimates
#' @param idx Indices to carry out estimation over
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param hal_ind A logical value indicating whether to estimate via HAL
#' @param sl_learners A character vector containing the names of the superlearner algorithms
#'
#' @return A list containing the estimate of E[Y|A=a,X,W]
#' @export
est_g <- function(idx,A, X, R, kappa_hat,
                  hal_ind,
                  sl_learners) {

  if (hal_ind==TRUE) { # estimate via HAL
    g_hat<- hal9001::fit_hal(Y=A[idx],X=X[idx,,drop=FALSE],weights=R[idx]/kappa_hat[idx])
  } else { # estimate via SL
    g_hat <- SuperLearner::SuperLearner(Y=A[idx],X=X[idx,,drop=FALSE],
                                        family=binomial(),SL.library=sl_learners,
                                        obsWeights=R[idx]/kappa_hat[idx])
  }
  return(g_hat)
}

# note: handle user supplied r probs outside of this function
est_kappa <- function (idx,Z, R,
                       hal_ind,
                       sl_learners) {

  # if estimating via highly adaptive lasso
  if (hal_ind==TRUE) { # estimate via HAL
    kappa_hat <- hal9001::fit_hal(Y=R[idx],X=Z[idx,], family = 'binomial')
    kappa_hat <- predict(kappa_hat, new_data=Z)
  }

  # if estimating via superlearner
  if (hal_ind==FALSE) { # estimate via SL
    loadNamespace("SuperLearner")
    kappa_hat <- SuperLearner::SuperLearner(Y=R[idx],X=Z[idx,],
                                            family=binomial(),SL.library=sl_learners)
    kappa_hat <- predict(kappa_hat, newdata=Z)$pred

  }

  return(kappa_hat)

}


#' @title est_varphi
#'
#' @description Function for obtaining estimate of E[phi|Z]
#' @param idx Indices to carry out estimation over
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param Z Dataframe containing all non-missing variables
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param phi_hat A list containing the estimate of phi
#' @param eem_ind A logical value indicating whether to estimate via EEM
#'
#' @return A list containing the estimate of E[phi|Z]
#' @export
est_varphi_main <- function(idx, R,Z,
                            phi_1_hat, phi_0_hat,
                            kappa_hat,
                            eem_ind,
                            hal_ind,
                            sl_learners){


  if (eem_ind==TRUE) { # estimate via EEM
    return(est_varphi_eem(idx, R, Z,
                          phi_1_hat, phi_0_hat,
                          kappa_hat,
                          hal_ind,
                          sl_learners))
  } else {
    # browser()
    return(est_varphi(idx, R, Z,
                      phi_1_hat, phi_0_hat,
                      hal_ind,
                      sl_learners))
  }


}

#' @title Estimation of E[phi|Z]
#' @description Outer function for obtaining estimate of E[phi|Z]
#' @param idx Indices to carry out estimation over
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param Z Dataframe containing all non-missing variables
#' @param phi_1_hat A numeric vector containing the fitted values of phi under A=1
#' @param phi_0_hat A numeric vector containing the fitted values of phi under A=0
#' @param hal_ind A logical value indicating whether to estimate via HAL
#' @param sl_learners A character vector containing the names of the superlearner algorithms
#' @return A list containing the estimate of E[phi|Z]
#' @export
est_varphi <- function(idx, R, Z,
                       phi_1_hat, phi_0_hat,
                       hal_ind,
                       sl_learners) {
  print('of is g')
  if (hal_ind==TRUE) { # estimate via HAL
    varphi_1_hat <- hal9001::fit_hal(Y=phi_1_hat[idx],X=Z[idx,,drop=FALSE],weights=R[idx])
    varphi_0_hat <- hal9001::fit_hal(Y=phi_0_hat[idx],X=Z[idx,,drop=FALSE],weights=R[idx])
  } else { # estimate via SL
    varphi_1_hat <- SuperLearner::SuperLearner(Y=phi_1_hat[idx],X=Z[idx,,drop=FALSE],
                                               family=gaussian(),SL.library=sl_learners,
                                               obsWeights=R[idx])
    varphi_0_hat <- SuperLearner::SuperLearner(Y=phi_0_hat[idx],X=Z[idx,,drop=FALSE],
                                               family=gaussian(),SL.library=sl_learners,
                                               obsWeights=R[idx])
  }

  return(list(varphi_1_hat=varphi_1_hat,varphi_0_hat=varphi_0_hat))
}


#' @title est_varphi_eem
#' @description Function for obtaining estimate of E[phi|Z] via empirical efficiency maximization
#'
#' @param data A data frame
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param Z Dataframe containing all non-missing variables
#' @param phi_1_hat Vector of predicted phi under A=1
#' @param phi_0_hat Vector of predicted phi under A=0
#' @param eem_ind A logical value indicating whether to estimate via EEM
#' @param hal_ind A logical value indicating whether to estimate via HAL
#' @param sl_learners A character vector containing the names of the superlearner algorithms
#'
#' @return A list containing the estimate of E[phi|Z]
#' @export
est_varphi_eem <- function(idx, R, Z,
                           phi_1_hat, phi_0_hat,
                           kappa_hat,
                           hal_ind,
                           sl_learners) {
  print('of is ye')
  # Make pseudo outcomes
  ytilde1 <- (R/kappa_hat -1)^(-1) * (R/kappa_hat) * phi_1_hat
  ytilde0 <- (R/kappa_hat -1)^(-1) * (R/kappa_hat) * phi_0_hat

  # Estimate E[phi|Z] via EEM
  if (hal_ind==TRUE) { # estimate via HAL
    varphi_1_hat <- hal9001::fit_hal(Y=ytilde1[idx],X=Z[idx,,drop=FALSE],
                                   family=gaussian(),SL.library=sl_learners,
                                   obsWeights=(R[idx]/kappa_hat[idx] - 1)^2)
    varphi_0_hat <- hal9001::fit_hal(Y=ytilde0[idx],X=Z[idx,,drop=FALSE],
                                   family=gaussian(),SL.library=sl_learners,
                                   obsWeights=(R[idx]/kappa_hat[idx] - 1)^2)

  } else { # estimate via SL
    varphi_1_hat <- SuperLearner::SuperLearner(Y=ytilde1[idx],X=Z[idx,,drop=FALSE],
                                             family=gaussian(),SL.library=sl_learners,
                                             obsWeights=(R[idx]/kappa_hat[idx] - 1)^2)
    varphi_0_hat <- SuperLearner::SuperLearner(Y=ytilde0[idx],X=Z[idx,,drop=FALSE],
                                               family=gaussian(),SL.library=sl_learners,
                                               obsWeights=(R[idx]/kappa_hat[idx] - 1)^2)

  }
  return(list(varphi_1_hat=varphi_1_hat,varphi_0_hat=varphi_0_hat))
}

