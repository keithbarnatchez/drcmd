#' @title Doubly-robust causal inference with missing data
#'
#'
#' @description Doubly-robust estimation of counterfactual means in the presence of missing
#  'data. The drcmd() function estimates counterfactual means for binary point treatments
#'  and reports average treatment effects, as well as causal risk ratios and odds ratios
#'  for binary outcomes. General missingness patterns in the data are allowed and automatically
#'  determined by the function -- the only requirement is that any missingness occurs at random
#'  conditional on variables that are always available. For scenarios where non-missingness
#'  probabilities are known, as is common in two-phase sampling designs, users can provide
#'  the non-missingness probabilities through the Rprobs argument. Users can fit nuisance
#'  functions through either highly-adaptive LASSO (HAL) or SuperLearner, the latter of which
#'  the user must specify libraries.
#'
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param W Dataframe containing variables solely predictive of missingness, but not
#' a cause of the outcome or exposure.
#' @param R (optional) A character string specifying the missingness indicator,
#' where 0 indicates missing data. If not specified, the function will create the
#' missingness indicator by identifying the missingness pattern in the data
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl_learners A character vector containing SuperLearner libraries to use
#' for estimating all nuisance functions. User can alternatively specify libraries
#' for each nuisance function for added flexibility
#' @param m_sl_learners A character vector containing SuperLearner libraries to
#' be used for the outcome regression
#' @param g_sl_learners A character vector containing SuperLearner libraries to be
#' used for propensity score estimation
#' @param r_sl_learners A character vector containing SuperLearner libraries to be
#' used for missingness indicator estimation
#' @param po_sl_learners A character string containing SuperLearner libraries to be
#' used for the pseudo-outcome vector
#' @param eem_ind A logical indicating whether to use ensemble of estimators
#' @param Rprobs A vector of probabilities for the missingness indicator. Only
#' suitable for study designs where researcher controls mechanism by which variables
#' are missing (e.g. two-phase sample designs). Defaults to NA, in which case
#' missingness probabilities are estimated.
#' @param k A numeric indicating the number of folds for cross-fitting
#' @param nboot A numeric indicating the number of desired bootstrap samples.
#' If >0, uses bootstrap to obtain SEs. If =0, uses asymptotic analytical SEs.
#' @return A list of results
#'
#' @examples
#' n <- 1e3
#' X <- rnorm(n) ; A <- rbinom(n,1,plogis(X)) ; Y <- rnorm(n) + A + X
#' Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) ; X <- as.data.frame(X)
#' Y[R==0] <- NA
#'
#' results <- drcmd(Y,A,X,R=R)
#'
#' @export
drcmd <- function(Y, A, X, W=NA, R=NA,
                  hal_ind=TRUE,
                  sl_learners=NULL,
                  m_sl_learners=NULL,g_sl_learners=NULL,r_sl_learners=NULL,po_sl_learners=NULL,
                  eem_ind=FALSE, Rprobs=NA, k=1,
                  nboot=0) {

  require(SuperLearner)

  # Throw errors if anything is entered incorrectly
  # check_entry_errors(Y,A,X,W,R, hal_ind,sl_learners,eem_ind,Rprobs)
  if (is.na(W)) {
    W <- X[,0]
  }

  # Identify missing data structure
  V <- find_missing_pattern(Y,A,X,W)
  Z <- V$Z ; R <- V$R ; X <- V$X ; Y <- V$Y ; A <- V$A

  # Estimate
  res <- list()
  res$results <- drcmd_est(Y,A,X,Z,R,hal_ind,sl_learners,eem_ind,Rprobs,k)

  # Additional packaging
  params <- list()
  params$hal_ind <- hal_ind
  params$sl_learners <- sl_learners
  params$m_sl_learners <- m_sl_learners
  params$g_sl_learners <- g_sl_learners
  params$r_sl_learners <- r_sl_learners
  params$po_sl_learners <- po_sl_learners
  params$eem_ind <- eem_ind
  params$Rprobs <- Rprobs
  params$k <- k
  res$Z <- colnames(Z)
  res$U <- V$U
  class(res) <- 'drcmd'

  # Return results
  return(res)

}

#' @title Obtain doubly-robust counterfactual mean estimates through cross-fitting
#'
#' @description Outer function for estimating counterfactual means through cross-fitting.
#' Calls \code{drcmd_est_fold} to obtain estimates for each cross-fitting fold
#' @param data A data frame
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param Z Dataframe containing all variables that are never subject to missingness
#' @param R Binary missingness indicator
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl_learners A character string indicating the superlearner library to use
#' @param eem_ind A logical indicating whether to use ensemble of estimators
#' @param Rprobs A vector of probabilities for R
#' @param k A numeric indicating the number of folds for cross-fitting
#'
#' @return A list of results
#' @export
drcmd_est <- function(Y,A,X,Z,R,hal_ind,sl_learners,eem_ind,Rprobs,k) {

  # Divide the data into k random "train-test" splits
  if (k>1) {
    splits <- create_folds(length(Y),k)

    # Use lapply to get point ests for each fold
    res <- lapply(1:k, function(i) drcmd_est_fold(splits[[i]],Y,A,X,Z,R,hal_ind,sl_learners,eem_ind,Rprobs))

    # Extract ests and SEs from each fold
    ests_df <- do.call(rbind, lapply(res, function(x) x$ests))
    vars_df <- do.call(rbind, lapply(res, function(x) x$vars))

    # Combine
    ests <- colMeans(ests_df)
    e_ests_sq <- colMeans( (sweep(ests_df,2,ests,'-'))^2 )
    vars <- colMeans( sweep(vars_df,2,e_ests_sq,'+') )
    res <- list(estimates=ests,ses=sqrt(vars))
    return(res)
  } else {
    splits <- create_folds(length(Y),k)
    res <- drcmd_est_fold(splits,Y,A,X,Z,R,hal_ind,sl_learners,eem_ind,Rprobs)
    res <- list(estimates=res$ests,ses=sqrt(res$vars))
    return(res)
  }
}

#' @title Calculate point estimates within a single cross-fitting fold
#'
#' @description Outer function for obtaining point estimates for single fold
#' @param splits A list of train/test splits
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param Z Dataframe containing all variables that are never subject to missingness
#' @param R Binary missingness indicator
#' @param hal_ind A logical indicating whether to use highly adaptive lasso
#' @param sl_learners A character string indicating the superlearner library to use
#' @param eem_ind A logical indicating whether to use ensemble of estimators
#' @param Rprobs A vector of probabilities for R
#' @param k A numeric indicating the number of folds for cross-fitting
#'
#' @return A list of results
#' @export
drcmd_est_fold <- function(splits,Y,A,X,Z,R,hal_ind,sl_learners,eem_ind,Rprobs) {

  # Get the training and test data
  train <- splits$train
  test <- splits$test

  # Get the nuisance estimates
  nuisance_ests <- get_nuisance_ests(train,Y,A,X,Z,R,hal_ind,sl_learners,Rprobs)

  # Form full data EIF, phi
  phi_hat <- get_phi_hat(Y,A,X,R,nuisance_ests$g_hat,nuisance_ests$m_a_hat,
                      nuisance_ests$kappa_hat,hal_ind,sl_learners)

  # Get varphi
  phi_1_hat <- phi_hat$phi_1_hat ; phi_0_hat <- phi_hat$phi_0_hat
  varphi_hat <- est_varphi_main(test,R,Z,phi_1_hat,phi_0_hat,
                                nuisance_ests$kappa_hat,
                                eem_ind,
                                hal_ind,sl_learners)

  # Form final est for this fold
  ests <- est_psi(test, R, Z, nuisance_ests$kappa_hat, phi_hat,varphi_hat)

  # Return results
  return(ests)

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
#' @param sl_learners A character string indicating the superlearner library to use
#' @return A list containing the estimated EIFs for A=1 and A=0
#'
#' @export
#'
get_phi_hat <- function(Y, A, X, R, g_hat, m_a_hat, kappa_hat,
                        hal_ind,sl_learners) {


  data <- cbind(X,A)

  # get predicted values for current data
  if (hal_ind==TRUE) { # estimate via HAL
    # get fitted outcome values under A=1 and A=0
    data[,ncol(data)] <- 1 ; m_1_hat <- predict(m_a_hat, new_data=data)
    data[,ncol(data)] <- 0 ; m_0_hat <- predict(m_a_hat, new_data=data)

    # get propensity scores
    g_hat <- predict(g_hat, new_data=data)
  } else { # estimate via SL *** need dynamic family ***
    data[,ncol(data)] <- 1 ; m_1_hat <- predict(m_a_hat, newdata=data)$pred
    data[,ncol(data)] <- 0 ; m_0_hat <- predict(m_a_hat, newdata=data)$pred

    # get propensity scores
    g_hat <- predict(g_hat, newdata=data)$pred
  }
  data[,ncol(data)] <- A

  # get fitted values under A=1
  phi_1_hat <- m_1_hat + A*(Y - m_1_hat)/g_hat
  phi_0_hat <- m_0_hat + (1-A)*(Y - m_0_hat)/(1-g_hat)

  # Set values where R==0 to 0
  phi_1_hat[R==0] <- 0
  phi_0_hat[R==0] <- 0

  return(list(phi_1_hat=phi_1_hat,phi_0_hat=phi_0_hat))
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
est_psi <- function(idx, R, Z,
                    kappa_hat, phi_hat,varphi_hat) {

  n <- length(idx)

  # Set up necessary objects
  phi_1_hat <- phi_hat$phi_1_hat
  phi_0_hat <- phi_hat$phi_0_hat
  varphi_1_hat <- predict(varphi_hat$varphi_1_hat, newdata=Z)$pred
  varphi_0_hat <- predict(varphi_hat$varphi_0_hat, newdata=Z)$pred

  # Form estimates of psi1 and psi0 via EICs
  psi_1_ic <- (R[idx]/kappa_hat[idx])*phi_1_hat[idx] - (R[idx]/kappa_hat[idx] - 1)*varphi_1_hat[idx]
  psi_0_ic <- (R[idx]/kappa_hat[idx])*phi_0_hat[idx] - (R[idx]/kappa_hat[idx] - 1)*varphi_0_hat[idx]

  psi_1_hat <- mean(psi_1_ic)
  psi_0_hat <- mean(psi_0_ic)

  # Get specific
  psi_hat_ate <- psi_1_hat - psi_0_hat
  psi_hat_rr <- psi_1_hat/psi_0_hat
  psi_hat_or <- (psi_1_hat/(1-psi_1_hat)) / (psi_0_hat/(1-psi_0_hat))

  return(list(ests = data.frame(psi_1_hat=mean(psi_1_ic),
                                psi_0_hat=mean(psi_0_ic),
                                psi_hat_ate=mean(psi_1_ic - psi_0_ic),
                                psi_hat_rr=mean(psi_1_ic)/mean(psi_0_ic),
                                psi_hat_or=(mean(psi_1_ic)/(1-mean(psi_1_ic))) / (mean(psi_0_ic)/(1-mean(psi_0_ic)))
              ),
              vars = data.frame(psi_1_hat=var(psi_1_ic)/n,
                                psi_0_hat=var(psi_0_ic)/n,
                                psi_hat_ate=var(psi_1_ic - psi_0_ic)/n,
                                psi_hat_rr=var(psi_1_ic)/var(psi_0_ic)/n,
                                psi_hat_or=((var(psi_1_ic)/(1-var(psi_1_ic))) / (var(psi_0_ic)/(1-var(psi_0_ic))))/n
              )
         )
  )

  # return(list(psi_1_hat=list(est=mean(psi_1_ic), var=var(psi_1_ic)),
  #             psi_0_hat=list(est=mean(psi_0_ic), var=var(psi_0_ic)),
  #             psi_hat_ate=list(est=mean(psi_1_ic - psi_0_ic), var=var(psi_1_ic - psi_0_ic))
  #       )
  # )

}
