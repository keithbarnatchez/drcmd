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
  if (any(is.na(Rprobs) )) {
    kappa_hat <- est_kappa(idx,Z,R,r_learners)
  }
  if (!is.null(cutoff)) {
    kappa_hat <- truncate_r(kappa_hat,cutoff)
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
  m_a_hat <- SuperLearner::SuperLearner(Y=Y[idx],
                                        X=data[idx,],
                                        SL.library=m_learners,
                                        family=yfam,
                                        # method=method,
                                        obsWeights=R[idx]/kappa_hat[idx])
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

  g_hat <- SuperLearner::SuperLearner(Y=A[idx],
                                      X=X[idx,,drop=FALSE],
                                      family=binomial(),
                                      SL.library=g_learners,
                                      # method='method.NNloglik',
                                      obsWeights=R[idx]/kappa_hat[idx])
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

  # browser()

  if (all(R==1)) { # when no missing data at all, no need for pseudo-outcome reg
    return(list(varphi_1_hat=0,
                varphi_0_hat=0,
                varphi_diff_hat=0))
  }

  if (eem_ind==TRUE) { # estimate via EEM
    return(est_varphi_eem(idx, R, Z,
                          phi_1_hat, phi_0_hat,
                          kappa_hat,
                          po_learners,
                          Y))
  } else {
    # browser()
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

  varphi_diff_hat <- SuperLearner::SuperLearner(Y=phi_1_hat[idx]-phi_0_hat[idx],X=Z[idx,,drop=FALSE],
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
  if (all(po_learners=='hal')) { # estimate via HAL
    varphi_1_hat <- hal9001::fit_hal(Y=ytilde1[idx],X=Z[idx,,drop=FALSE],
                                   weights=(R[idx]/kappa_hat[idx] - 1)^2)
    varphi_0_hat <- hal9001::fit_hal(Y=ytilde0[idx],X=Z[idx,,drop=FALSE],
                                   weights=(R[idx]/kappa_hat[idx] - 1)^2)

    varphi_diff_hat <- hal9001::fit_hal(Y=ytilde1[idx]-ytilde0[idx],X=Z[idx,,drop=FALSE],
                                     weights=(R[idx]/kappa_hat[idx] - 1)^2)

    varphi_1_hat <- predict(varphi_1_hat, new_data=Z)
    varphi_0_hat <- predict(varphi_0_hat, new_data=Z)
    varphi_diff_hat <- predict(varphi_diff_hat, new_data=Z)

  } else { # estimate via SL
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
  }
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


SL.earth_weight <- function(Y, X, newX, family, obsWeights, id,
                            degree = 2, penalty = 3,
                            nk = max(21, 2 * ncol(X) + 1),
                            pmethod = "backward", nfold = 0,
                            ncross = 1, minspan = 0, endspan = 0, ...) {
  require("earth")

  # Fit model with weights
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(
      x = X, y = Y, weights = obsWeights,
      degree = degree, nk = nk, penalty = penalty,
      pmethod = pmethod, nfold = nfold, ncross = ncross,
      minspan = minspan, endspan = endspan
    )
  } else if (family$family == "binomial") {
    fit.earth <- earth::earth(
      x = X, y = Y, weights = obsWeights,
      degree = degree, nk = nk, penalty = penalty,
      pmethod = pmethod, nfold = nfold, ncross = ncross,
      minspan = minspan, endspan = endspan,
      glm = list(family = binomial)
    )
  } else {
    stop("SL.earth_weight currently supports only gaussian and binomial families.")
  }

  # Predictions
  pred <- predict(fit.earth, newdata = newX, type = "response")

  # Return model object
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.earth_weight"

  return(out)
}

SL.gam_interact <- function(Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4, ...) {
  requireNamespace("mgcv", quietly = TRUE)

  # Identify continuous vs categorical variables
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))

  # Construct formula with smooth main effects and tensor interactions
  if (sum(cts.x) > 1) {  # Only include interactions if >1 continuous variable
    gam.model <- as.formula(paste(
      "Y ~",
      paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ", k=5, bs='tp')", sep = ""), collapse = " + "),
      "+ te(", paste(colnames(X[, cts.x, drop = FALSE]), collapse = ", "), ", bs='tp')"
    ))
  } else if (sum(!cts.x) > 0) {  # If only categorical variables
    gam.model <- as.formula(paste("Y ~", paste(colnames(X), collapse = "+")))
  } else {  # If only one continuous variable, no interactions possible
    gam.model <- as.formula(paste(
      "Y ~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ", k=5, bs='tp')", sep = ""), collapse = " + ")
    ))
  }

  # Fit GAM model using mgcv (supports `te()`)
  fit.gam <- mgcv::gam(
    gam.model, data = data.frame(Y = Y, X), family = family,
    weights = obsWeights, method = "REML"
  )

  # Predict on new data
  pred <- predict(fit.gam, newdata = newX, type = "response")

  # Return model in required format
  fit_obj <- list(object = fit.gam)
  class(fit_obj) <- "SL.gam_interact"

  return(list(pred = pred, fit = fit_obj))
}

# **Prediction wrapper for SuperLearner**
predict.SL.gam_interact <- function(object, newdata, ...) {
  predict(object$object, newdata = newdata, type = "response")
}




