#' @title Obtain nuisance function estimates
#'
#' @description Function for obtaining estimates of all relevant nuisance functions
#' @param idx Indices to carry out estimation over
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param Z Dataframe containing all non-missing variables
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param m_learners A character vector containing learners to be used for the
#' outcome regression. User can specify 'hal' or a vector of SuperLearner libraries
#' @param g_learners A character vector containing learners to be used for the
#' propensity score. User can specify 'hal' or a vector of SuperLearner libraries
#' @param r_learners A character vector containing learners to be used for the
#' missingness indicator regression. User can specify 'hal' or a vector of
#' SuperLearner libraries
#' @param eem_ind A logical indicating whether to use empirical efficiency maximization
#' @param Rprobs A vector of probabilities for R
#'
#' @return A list of nuisance estimates
#' @export
get_nuisance_ests <- function(idx,Y,A,X,Z,R,
                              m_learners,g_learners,r_learners,
                              Rprobs,cutoff) {

  kappa_hat <- Rprobs
  if (all(R==1)) { # if no missing data
    kappa_hat <- rep(1,nrow(X))
  } else {
    if (any(is.na(Rprobs) )) {
      kappa_hat <- est_kappa(idx,Z,R,r_learners)
    }
    if (!is.null(cutoff)) {
      kappa_hat <- truncate_r(kappa_hat,cutoff)
    }
  }
  m_a_hat <- est_m_a(idx,Y,A,X,R,kappa_hat,m_learners)
  g_hat <- est_g(idx,A,X,R,kappa_hat,g_learners)

  # Trimming
  if (!is.null(cutoff)) {
    g_hat <- truncate_g(g_hat,cutoff)
  }

  return(list(kappa_hat=kappa_hat,m_a_hat=m_a_hat,g_hat=g_hat))

}


#' @title Estimating outcome regression
#'
#' @description Function for obtaining estimate of E[Y|A=a,X]. Estimation is
#' carried out with SuperLearer, using learners specified in m_learners
#' @param idx Indices to carry out estimation over
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param kappa_hat A numeric vector containing the fitted values of kappa
#' @param m_learners A character vector containing the names of the superlearner algorithms
#' @return A list containing the estimate of E[Y|A=a,X]
#' @export
#' @examples
#' \dontrun{
#' n <- 1000
#' X <- rnorm(n)
#' A <- rbinom(n,1,plogis(X))
#' R <- rbinom(n,1,plogis(X))
#' X <- data.frame(X)
#' Y <- A + X rnorm(n)
#' m_learners <- c('SL.glm','SL.gam')
#' est_m_a(idx=1:n, A=A, X=X, R=R, kappa_hat=kappa_hat, g_learners=g_learners)
#' }
est_m_a <- function(idx, Y, A, X, R,
                    kappa_hat,
                    m_learners) {

  yfam <- gaussian()
  # method <- 'method.NNLS'
  if (check_binary(Y)) { # will treat [0,1] bounded outcome as binary
    yfam <- binomial()
    # method <- 'method.NNloglik'
  }

  data <- cbind(X,A)

  # if outcome is binary, weights will generate unecessary warnings
  m_a_hat <- g_hat <- withCallingHandlers(
    SuperLearner::SuperLearner(Y=Y[idx],
                              X=data[idx,],
                              SL.library=m_learners,
                              family=yfam,
                              # method=method,
                              obsWeights=R[idx]/kappa_hat[idx]),
    warning = function(w) {
      if (conditionMessage(w) == "non-integer #successes in a binomial glm!") {
        invokeRestart("muffleWarning")
      }
    })

  # get fitted values under A=1
  data[,ncol(data)] <- 1
  m_1_hat <- predict(m_a_hat, newdata=data)$pred
  data[,ncol(data)] <- 0
  m_0_hat <- predict(m_a_hat, newdata=data)$pred

  return(list(m_1_hat=m_1_hat,m_0_hat=m_0_hat))
}


#' @title Estimate the propensity score
#'
#' @description Function for obtaining propensity score estimates P(A=a|X)
#' @param idx Indices to carry out estimation over
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param g_learners A character vector containing the names of the learners for estimation
#'
#' @return A list containing the estimate of E[Y|A=a,X,W]
#' @export
#' @examples
#' \dontrun{
#' n <- 1000
#' X <- rnorm(n)
#' A <- rbinom(n,1,plogis(X))
#' R <- rbinom(n,1,plogis(X))
#' X <- data.frame(X)
#' g_learners <- c('SL.glm','SL.gam')
#' est_g(idx=1:n, A=A, X=X, R=R, kappa_hat=kappa_hat, g_learners=g_learners)
#' }
est_g <- function(idx,A, X, R, kappa_hat,
                  g_learners) {

  # regression with weights + binary outcome can create unnecessary warnings
  g_hat <- withCallingHandlers(
    SuperLearner::SuperLearner(
      Y = A[idx],
      X = X[idx, , drop = FALSE],
      family = binomial(),
      SL.library = g_learners,
      obsWeights = R[idx] / kappa_hat[idx]
    ),
    warning = function(w) {
      if (conditionMessage(w) == "non-integer #successes in a binomial glm!") {
        invokeRestart("muffleWarning")
      }
    }
  )

  g_hat <- predict(g_hat, newdata=X)$pred

  return(g_hat)
}

#' @title Estimate complete case propensity scores
#'
#' @description Function for obtaining complete case propensity score estimates
#'
#' @param idx Indices to carry out estimation over
#' @param Z Dataframe containing all non-missing variables
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param r_learners A character vector specifying learners to be used for estimation
#'
#' @return A numeric vector containing the estimate of E[R|Z]
#' @export
#' @examples
#' \dontrun{
#' n <- 1000
#' Z <- rnorm(n)
#' R <- rbinom(n,1,plogis(Z))
#' Z <- data.frame(Z)
#' g_learners <- c('SL.glm','SL.gam')
#' est_g(idx=1:n, Z=Z, R=R, r_learners=r_learners)
#' }
est_kappa <- function (idx,Z, R,
                       r_learners) {

  # if estimating via highly adaptive lasso
  if (all(r_learners=='hal')) { # estimate via HAL
    kappa_hat <- hal9001::fit_hal(Y=R[idx],X=Z[idx,], family = 'binomial')
    kappa_hat <- predict(kappa_hat, new_data=Z)
  } else {  # estimate via SL
    loadNamespace("SuperLearner")
    kappa_hat <- SuperLearner::SuperLearner(Y=R[idx],X=Z[idx,,drop=FALSE],
                                            family=binomial(),
                                            SL.library=r_learners)
    kappa_hat <- predict(kappa_hat, newdata=Z)$pred
  }

  return(kappa_hat)

}

#' @title Perform pseudo-outcome regression
#'
#' @description Function for obtaining estimate of E[phi_a|Z] through pseudo-outcome
#' regression. Calls inner function to perform estimation, depending on whether
#' user wishes to perform empirical efficiency maximization or not
#' @param idx Indices to carry out estimation over
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param Z Dataframe containing all non-missing variables
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param phi_hat A list containing the estimate of phi
#' @param eem_ind A logical value indicating whether to estimate via EEM
#'
#' @return A list containing the estimate of E[phi_a|Z] for a=0 and a=1
#' @export
est_varphi_main <- function(idx, R,Z,
                            phi_1_hat, phi_0_hat,
                            kappa_hat,
                            eem_ind,
                            po_learners,
                            Y){

  if (all(R==1)) { # when no missing data at all, no need for pseudo-outcome reg
    return(list(varphi_1_hat=rep(0,length(kappa_hat)),
                varphi_0_hat=rep(0,length(kappa_hat)),
                varphi_diff_hat=rep(0,length(kappa_hat))))
  }

  if (eem_ind==TRUE) { # estimate via EEM
    return(est_varphi_eem(idx, R, Z,
                          phi_1_hat, phi_0_hat,
                          kappa_hat,
                          po_learners,
                          Y))
  } else {
    return(
      est_varphi(idx, R, Z,
                      phi_1_hat, phi_0_hat,
                      po_learners,
                      Y)
           )
  }


}

#' @title Perform pseudo-outcome regression with conventional loss function
#' @description Outer function for obtaining estimate of E[phi|Z]
#' @param idx Indices to carry out estimation over
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param Z Dataframe containing all non-missing variables
#' @param phi_1_hat A numeric vector containing the fitted values of phi under A=1
#' @param phi_0_hat A numeric vector containing the fitted values of phi under A=0
#' @param po_learners A character vector containing the names of the superlearner algorithms
#' @return A list containing the estimate of E[phi|Z] for a=0 and a=1
#' @export
est_varphi <- function(idx, R, Z,
                       phi_1_hat, phi_0_hat,
                       po_learners,
                       Y) {

  varphi_diff_hat <- SuperLearner::SuperLearner(Y=phi_1_hat[idx]-phi_0_hat[idx],
                                                X=Z[idx,,drop=FALSE],
                                                family=gaussian(),
                                                SL.library=po_learners,
                                                obsWeights=R[idx])

  # if Y binary or all in (0,1), shift and scale IFs
  fam <- gaussian()
  if (check_binary(Y)) { # will treat [0,1] bounded outcome as binary
    fam <- binomial()

    max1 <- max(phi_1_hat)
    min1 <- min(phi_1_hat)
    max0 <- max(phi_0_hat)
    min0 <- min(phi_0_hat)
    mindiff <- min(phi_1_hat - phi_0_hat)
    maxdiff <- max(phi_1_hat - phi_0_hat)

    phi_1_hat <- (phi_1_hat - min1)/(max1 - min1)
    phi_0_hat <- (phi_0_hat - min0)/(max0 - min0)
    phi_diff_hat <- (phi_1_hat - phi_0_hat - mindiff)/(maxdiff - mindiff)

  }

  varphi_1_hat <- SuperLearner::SuperLearner(Y=phi_1_hat[idx],X=Z[idx,,drop=FALSE],
                                             family=fam,
                                             SL.library=po_learners,
                                             obsWeights=R[idx])
   varphi_0_hat <- SuperLearner::SuperLearner(Y=phi_0_hat[idx],X=Z[idx,,drop=FALSE],
                                             family=fam,
                                             SL.library=po_learners,
                                             obsWeights=R[idx])

  varphi_1_hat <- predict(varphi_1_hat, newdata=Z)$pred
  varphi_0_hat <- predict(varphi_0_hat, newdata=Z)$pred
  varphi_diff_hat <- predict(varphi_diff_hat, newdata=Z)$pred

  # rescale varphi_1_hat and varphi_0_hat if Y was binary
  if (check_binary(Y)) {
    varphi_1_hat <- varphi_1_hat * (max1 - min1) + min1
    varphi_0_hat <- varphi_0_hat * (max0 - min0) + min0
  }


  return(list(varphi_1_hat=varphi_1_hat,varphi_0_hat=varphi_0_hat,
              varphi_diff_hat=varphi_diff_hat))
}


#' @title Perform pseudo-outcome regression with empirical efficiency maximization
#' @description Function for obtaining estimate of E[phi_a|Z] via empirical efficiency
#' maximization
#'
#' @param data A data frame
#' @param R Binary missingness indicator, where 0 indicates missing data
#' @param Z Dataframe containing all non-missing variables
#' @param phi_1_hat Vector of predicted phi under A=1
#' @param phi_0_hat Vector of predicted phi under A=0
#' @param eem_ind A logical value indicating whether to estimate via EEM
#' @param po_learners A character vector containing the names of the superlearner algorithms
#'
#' @return A list containing the estimate of E[phi_a|Z] for a=0 and a=1
#' @export
est_varphi_eem <- function(idx, R, Z,
                           phi_1_hat, phi_0_hat,
                           kappa_hat,
                           po_learners,
                           Y) {

  # Make pseudo outcomes
  ytilde1 <- (R/kappa_hat -1)^(-1) * (R/kappa_hat) * phi_1_hat
  ytilde0 <- (R/kappa_hat -1)^(-1) * (R/kappa_hat) * phi_0_hat

  # Estimate E[phi|Z] via EEM
  #### ***** need to update hal code below
  varphi_1_hat <- SuperLearner::SuperLearner(Y=ytilde1[idx],X=Z[idx,,drop=FALSE],
                                           family=gaussian(),SL.library=po_learners,
                                           obsWeights=(R[idx]/kappa_hat[idx] - 1)^2)
  varphi_0_hat <- SuperLearner::SuperLearner(Y=ytilde0[idx],X=Z[idx,,drop=FALSE],
                                             family=gaussian(),SL.library=po_learners,
                                             obsWeights=(R[idx]/kappa_hat[idx] - 1)^2)
  varphi_1_hat <- predict(varphi_1_hat, newdata=Z)$pred
  varphi_0_hat <- predict(varphi_0_hat, newdata=Z)$pred

  varphi_diff_hat <- SuperLearner::SuperLearner(Y=ytilde1[idx]-ytilde0[idx],X=Z[idx,,drop=FALSE],
                                                family=gaussian(),SL.library=po_learners,
                                                obsWeights=(R[idx]/kappa_hat[idx] - 1)^2)
  varphi_diff_hat <- predict(varphi_diff_hat, newdata=Z)$pred

  return(list(varphi_1_hat=varphi_1_hat,varphi_0_hat=varphi_0_hat,
              varphi_diff_hat=varphi_diff_hat))
}

#' @title SuperLearner wrapper for the highly-adaptive lasso
#'
#' @description Wrapper for the highly-adaptive lasso (HAL) algorithm
#'implemented in the hal9001 package
#' @param Y Outcome variable
#' @param X Matrix of covariates
#' @param newX Matrix of covariates for prediction
#' @param family Family of the outcome variable
#' @param obsWeights Observation weights
#' @param ... Additional arguments to pass to the HAL function
#' @return A list containing the prediction and the fitted model
#' @export
SL.hal9001 <- function(Y, X, newX, family, obsWeights, ...) {
  # Fit HAL model
  fit <- hal9001::fit_hal(X = X, Y = Y, family = family$family,
                          max_degree=2,num_knots=3,
                          reduce_basis=TRUE,
                          lambda=0.01,
                          weights = obsWeights, ...)

  # Predict on new data
  pred <- predict(fit, new_data = newX)

  # Return in required format
  fit_obj <- list(object = fit)
  class(fit_obj) <- "SL.hal9001"
  return(list(pred = pred, fit = fit_obj))
}

#' @title Prediction wrapper for the highly-adaptive lasso
#' @description Prediction wrapper for the highly-adaptive lasso (HAL) algorithm
#' implemented in the hal9001 package
#' @param Y Outcome variable
#' @param X Matrix of covariates
#' @param newX Matrix of covariates for prediction
#' @param family Family of the outcome variable
#' @param obsWeights Observation weights
#' @param ... Additional arguments to pass to the HAL function
#' @return Predicted values on new data
#' @export
predict.SL.hal9001 <- function(object, newdata, ...) {
  predict(object$object, new_data = newdata)
}


