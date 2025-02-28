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
#' @param W (optional) Dataframe containing variables solely predictive of missingness,
#'  but not a cause of the outcome or exposure.
#' @param R (optional) A character string specifying the missingness indicator,
#' where 0 indicates missing data. If not specified, the function will create the
#' missingness indicator by identifying the missingness pattern in the data
#' @param deafult_learners A character vector containing SuperLearner libraries to use
#' for estimating all nuisance functions. User can alternatively specify libraries
#' for each nuisance function for added flexibility
#' @param m_learners A character vector containing learners to be used for the
#' outcome regression. User can specify 'hal' or a vector of SuperLearner libraries
#' @param g_learners A character vector containing learners to be used for the
#' propensity score. User can specify 'hal' or a vector of SuperLearner libraries
#' @param r_learners A character vector containing learners to be used for the
#' missingness indicator regression. User can specify 'hal' or a vector of
#' SuperLearner libraries
#' @param po_learners A character vector containing learners to be used for the
#' pseudo-outcome regression. User can specify 'hal' or a vector of SuperLearner libraries
#' @param eem_ind A logical indicating whether to use empirical efficiency maximization
#' @param Rprobs A vector of probabilities for the missingness indicator. Only
#' suitable for study designs where researcher controls mechanism by which variables
#' are missing (e.g. two-phase sample designs). Defaults to NA, in which case
#' missingness probabilities are estimated.
#' @param k A numeric indicating the number of folds for cross-fitting
#' @param c Cutoff for treatment and complete case propensity scores. Estimates outside
#' of [c, 1-c] are set to c or 1-c, respectively
#' @param nboot A numeric indicating the number of desired bootstrap samples.
#' If >0, uses bootstrap to obtain SEs. If =0, uses asymptotic analytical SEs.
#' @return An S3 object of class drcmd containing estimation results, information
#' on the missing data structure, and parameters used in the estimation
#'
#' @examples
#' \dontrun{
#' n <- 2500
#' X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
#' Y <-  rbinom(n,1,plogis(X-A)) # rnorm(n) + A + X + X^2 + A*X + sin(X) # note: true ATE is 1
#' Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) # error-prone outcome measurements
#'
#' # Make A NA if R==0
#' A[R==0] <- NA
#' covariates <- data.frame(X1=X,X2=X2)
#'
#' # Obtain ATE estimates, fitting all nuisance models with ensemble of splines +
#' # GAMs (save for the pseudo-outcome regression, which is done with XGboost)
#' drcmd_res <- drcmd(Y,A,covariates,
#'                    default_learners= c('SL.gam','SL.earth'),
#'                    po_learners = 'SL.gam',
#'                    k=1,
#'                    eem_ind=F)
#' }
#' @export
drcmd <- function(Y, A, X, W=NA, R=NA,
                  default_learners=NULL,
                  m_learners=NULL,g_learners=NULL,r_learners=NULL,po_learners=NULL,
                  eem_ind=FALSE, Rprobs=NA, k=1, cutoff=0.01,
                  nboot=0) {

  require(SuperLearner)

  # Proxy variables may be empty
  if (any(is.na(W))) {
    W <- X[,0]
  }
  # browser()
  # Clean up learners
  learners <- clean_learners(default_learners,m_learners,g_learners,r_learners,
                             po_learners)

  # Throw errors if anything is entered incorrectly
  check_entry_errors(Y,A,X,W,R,eem_ind,Rprobs,k,nboot)

  # Identify missing data structure
  V <- find_missing_pattern(Y,A,X,W)
  Z <- V$Z ; R <- V$R ; X <- V$X ; Y <- V$Y ; A <- V$A

  # Check if Y is binary (all NAs set to 0 by this point)
  y_bin <- check_binary(Y)

  # Obtain estimates
  res <- list()
  res$results <- drcmd_est(Y,A,X,Z,R,
                           learners$m_learners,learners$g_learners,
                           learners$r_learners,learners$po_learners,
                           eem_ind,Rprobs,k,cutoff,y_bin)


  if('y' %in% colnames(Z)) {
    colnames(Z)[colnames(Z) == 'y'] <- 'Y'
  }

  # Additional packaging of params used in estimation
  params <- list()
  params$default_learners <- learners$default_learners
  params$m_learners <- learners$m_learners
  params$g_learners <- learners$g_learners
  params$r_learners <- learners$r_learners
  params$po_learners <- learners$po_learners
  params$eem_ind <- eem_ind
  params$Rprobs <- Rprobs
  params$k <- k
  res$Z <- colnames(Z)
  res$U <- V$U
  res$R <- R
  res$params <- params
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
#' @param m_learners A character vector containing learners to be used for the
#' outcome regression. User can specify 'hal' or a vector of SuperLearner libraries
#' @param g_learners A character vector containing learners to be used for the
#' propensity score. User can specify 'hal' or a vector of SuperLearner libraries
#' @param r_learners A character vector containing learners to be used for the
#' missingness indicator regression. User can specify 'hal' or a vector of
#' SuperLearner libraries
#' @param po_learners A character vector containing learners to be used for the
#' pseudo-outcome regression. User can specify 'hal' or a vector of SuperLearner libraries
#' @param eem_ind A logical indicating whether to use empirical efficiency maximization
#' @param Rprobs A vector of probabilities for R
#' @param k A numeric indicating the number of folds for cross-fitting
#'
#' @return A list of point estimates and standard errors for all estimands considered
#' @export
drcmd_est <- function(Y,A,X,Z,R,
                      m_learners,g_learners,r_learners,po_learners,
                      eem_ind,Rprobs,k,cutoff,y_bin) {

  # if (tml_ind) { # If estimating via TML
  #   1
  # }

  # Divide the data into k random "train-test" splits
  if (k>1) { # if doing cross-fitting
    splits <- create_folds(length(Y),k)

    # Use lapply to get results for each fold
    res <- lapply(1:k, function(i) drcmd_est_fold(splits[[i]],Y,A,X,Z,R,
                                                  m_learners,g_learners,r_learners,po_learners,
                                                  eem_ind,Rprobs,cutoff,y_bin))

    # Extract ests and SEs from each fold
    ests_df <- do.call(rbind, lapply(res, function(x) x$ests))
    vars_df <- do.call(rbind, lapply(res, function(x) x$vars))
    average_nuis <- as.data.frame(Reduce(`+`, lapply(res, `[[`, "nuis")) / length(res))

    # Combine estimates across folds
    ests <- colMeans(ests_df)
    if (FALSE) {
      e_ests_sq <- colMeans( (sweep(ests_df,2,ests,'-'))^2 )
      vars <- as.data.frame(t(colMeans( sweep(vars_df,2,e_ests_sq,'+') )))
    }
    else {
      vars <- est_ses_crossfit(res$ics)
    }
    # Get variances via fold-wise contributions to overall IC
    res <- list(estimates=as.data.frame(t(ests)),ses=sqrt(vars),nuis=average_nuis)

    return(res)
  } else { # if not doing cross-fitting
    splits <- create_folds(length(Y),k)
    res <- drcmd_est_fold(splits,Y,A,X,Z,R,
                          m_learners,g_learners,r_learners,po_learners,
                          eem_ind,Rprobs,cutoff,y_bin)
    res <- list(estimates=res$ests,ses=sqrt(res$vars),nuis=res$nuis)
    return(res)
  }
}

#' @title Cross-fit SE estimation via influence curve contributions
#'
#' @description Estimates cross-fit SEs by taking contributions to influence curve
#' of each fold
#'
#' @param ics A list of influence curves for each fold, obtained within drcmd_est
#'
#' @return A list of cross-fit SEs
#' @export
est_ses_crossfit <- function(ics) {

  # ics is a list of dataframes. rbind all the dataframes together
  ics_df <- do.call(rbind, ics)

  psi_1_ic <- ics_df$psi_1_ic
  psi_0_ic <- ics_df$psi_0_ic
  psi_ate_ic <- ics_df$psi_ate_ic

  psi_hat_rr_se <- psi_hat_or_se <- NA
  if (y_bin) { # if outcome is binary
    # Get estimates of E[Y(1)] and E[Y(0)] (needed for delta method SEs)
    psi_1_hat <- mean(psi_1_ic)
    psi_0_hat <- mean(psi_0_ic)

    psi_hat_rr <- psi_1_hat/psi_0_hat # risk ratio (only useful if binary)
    psi_hat_or <- (psi_1_hat/(1-psi_1_hat)) / (psi_0_hat/(1-psi_0_hat)) # odds ratio (only useful if binary)

    # standard errors (closed form via delta method)
    psi_hat_rr_se <- (Sig[1,1]/psi_0_hat^2 - 2*Sig[1,2]*psi_1_hat/(psi_0_hat^3) + Sig[2,2]*psi_1_hat^2/psi_0_hat^4)/n
    psi_hat_or_se <- psi_hat_or^2*(Sig[1,1]/(psi_1_hat^2*(1-psi_1_hat)^2) +
                                     Sig[2,2]/(psi_0_hat^2*(1-psi_0_hat)^2) -
                                     2*Sig[1,2]/(psi_1_hat*(1-psi_1_hat)*psi_0_hat*(1-psi_0_hat)))/n
  }


  return(
    data.frame(psi_1_hat=var(psi_1_ic)/n,
                    psi_0_hat=var(psi_0_ic)/n,
                    psi_hat_ate=var(psi_1_ic - psi_0_ic)/n,
                    psi_hat_ate_direct=var(psi_ate_ic)/n,
                    psi_hat_rr=psi_hat_rr_se,
                    psi_hat_or=psi_hat_or_se
  )
  )

}

#' @title Calculate point estimates within a single cross-fitting fold
#'
#' @description Outer function for obtaining point estimates for single fold
#' @param splits A list of train/test indices
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param Z Dataframe containing all variables that are never subject to missingness
#' @param R Binary missingness indicator
#' @param m_learners A character vector containing learners to be used for the
#' outcome regression. User can specify 'hal' or a vector of SuperLearner libraries
#' @param g_learners A character vector containing learners to be used for the
#' propensity score. User can specify 'hal' or a vector of SuperLearner libraries
#' @param r_learners A character vector containing learners to be used for the
#' missingness indicator regression. User can specify 'hal' or a vector of
#' SuperLearner libraries
#' @param po_learners A character vector containing learners to be used for the
#' pseudo-outcome regression. User can specify 'hal' or a vector of SuperLearner libraries
#' @param eem_ind A logical indicating whether to use empirical efficiency maximization
#' @param Rprobs A vector of probabilities for R
#' @param c
#'
#' @return A list of results from estimation on current fold
#' @export
drcmd_est_fold <- function(splits,Y,A,X,Z,R,
                           m_learners,g_learners,r_learners,po_learners,
                           eem_ind,Rprobs,cutoff,y_bin) {

  # Get training and test data indices
  train <- splits$train
  test <- splits$test

  # Get nuisance estimates
  nuisance_ests <- get_nuisance_ests(train,Y,A,X,Z,R,
                                     m_learners,g_learners,r_learners,
                                     Rprobs,cutoff)

  # Form full data EIF, phi
  phi_hat <- get_phi_hat(Y,A,X,R,
                      nuisance_ests$g_hat,nuisance_ests$m_a_hat,
                      nuisance_ests$kappa_hat)

  # Estimate varphi via pseudo-outcome regression
  phi_1_hat <- phi_hat$phi_1_hat
  phi_0_hat <- phi_hat$phi_0_hat
  varphi_hat <- est_varphi_main(test,R,Z,phi_1_hat,phi_0_hat,
                                nuisance_ests$kappa_hat,
                                eem_ind,
                                po_learners)

  # Form final estimate for this fold
  ests <- est_psi(test, R, Z, nuisance_ests$kappa_hat, phi_hat,varphi_hat,y_bin)

  # Additionally return nuisance est values
  nuis <- data.frame(kappa_hat=nuisance_ests$kappa_hat,
                     m_1_hat=nuisance_ests$m_a_hat$m_1_hat,
                     m_0_hat=nuisance_ests$m_a_hat$m_0_hat,
                     g_hat=nuisance_ests$g_hat,
                     phi_1_hat=phi_1_hat,
                     phi_0_hat=phi_0_hat,
                     varphi_1_hat=varphi_hat$varphi_1_hat,
                     varphi_0_hat=varphi_hat$varphi_0_hat)
  ests$nuis <- nuis

  # Return results
  return(ests)

}

#' @title Estimate causal effects through targeted maximum likelihood
#'
#' @description Function for obtaining estimates of causal estimands through
#' targeted maximum likelihood. To be called within drcmd_est_fold(). Updates
#' estimate of P(R=1|Z) and the plug-in estimator in a manner which removes the
#' plug-in bias incurred by an initial plug-in estimator
#'
#'
#'
#'
est_psi_tml <- function(idx, Y,A,X,
                        R, Z,
                        m_1_hat, m_0_hat, g_hat,
                        kappa_hat,varphi_hat) {

  # Set up necessary objects
  phi_1_hat <- phi_hat$phi_1_hat # full-data eif under A=1
  phi_0_hat <- phi_hat$phi_0_hat # full-data eif under A=0
  varphi_1_hat <- varphi_hat$varphi_1_hat # E[full-data eif under A=1|Z,R=1]
  varphi_0_hat <- varphi_hat$varphi_0_hat # E[full-data eif under A=0|Z,R=1]
  varphi_diff_hat <- varphi_hat$varphi_diff_hat

  # Update estimate of P(R=1|Z)





}

#' @title Contruct full data EIF estimate
#'
#' @description Function for constructing estimate of the full data EIF, given
#' propensity score and outcome model fits
#' @param Y Outcome variable. Can be continuous or binary
#' @param A A binary treatment variable (1=treated, 0=control)
#' @param X Dataframe containing baseline covariates
#' @param Z Dataframe containing all variables that are never subject to missingness
#' @param R Binary missingness indicator
#' @param g_hat Propensity score predictions
#' @param m_a_hat List of outcome predictions under A=0/1
#' @param kappa_hat Missingness probabilities
#'
#' @export
#'
get_phi_hat <- function(Y, A, X, R, g_hat, m_a_hat, kappa_hat) {

  m_1_hat <- m_a_hat$m_1_hat
  m_0_hat <- m_a_hat$m_0_hat

  # get fitted values under A=1 and A=0
  phi_1_hat <- m_1_hat + A*(Y - m_1_hat)/g_hat
  phi_0_hat <- m_0_hat + (1-A)*(Y - m_0_hat)/(1-g_hat)

  # Set values where R==0 to 0
  phi_1_hat[R==0] <- 0
  phi_0_hat[R==0] <- 0

  return(list(phi_1_hat=phi_1_hat,phi_0_hat=phi_0_hat))
}

#' @title Obtain estimates for current cross fitting fold
#' @description Function for obtaining counterfactual mean estimates for the current
#' fold of the cross-fitting procedure
#' @param idx A data frame
#' @param R A character string containing randomization variable name
#' @param Z A character vector containing the names of the variables in Z
#' @param kappa_hat A numeric vector containing the fitted values of kappa
#' @param phi_hat A list containing the estimate of E[phi|Z]
#' @param varphi_hat A list containing the estimate of E[phi|Z]
#'
#' @return A list of point estimates and standard errors for counterfactual means
#' and various counterfactual contrasts
#' @export
est_psi <- function(idx, R, Z,
                    kappa_hat, phi_hat,varphi_hat,y_bin) {

  # browser()
  n <- length(idx)

  # Set up necessary objects
  phi_1_hat <- phi_hat$phi_1_hat # full-data eif under A=1
  phi_0_hat <- phi_hat$phi_0_hat # full-data eif under A=0
  varphi_1_hat <- varphi_hat$varphi_1_hat # E[full-data eif under A=1|Z,R=1]
  varphi_0_hat <- varphi_hat$varphi_0_hat # E[full-data eif under A=0|Z,R=1]
  varphi_diff_hat <- varphi_hat$varphi_diff_hat

  # Form estimates of psi1 and psi0 via EICs
  psi_1_ic <- (R[idx]/kappa_hat[idx])*phi_1_hat[idx] - (R[idx]/kappa_hat[idx] - 1)*varphi_1_hat[idx]
  psi_0_ic <- (R[idx]/kappa_hat[idx])*phi_0_hat[idx] - (R[idx]/kappa_hat[idx] - 1)*varphi_0_hat[idx]
  psi_ate_ic <- (R[idx]/kappa_hat[idx])*(phi_1_hat[idx] - phi_0_hat[idx]) - (R[idx]/kappa_hat[idx] - 1)*varphi_diff_hat[idx]

  # Get estimates of E[Y(1)] and E[Y(0)]
  psi_1_hat <- mean(psi_1_ic)
  psi_0_hat <- mean(psi_0_ic)

  # Get estimates of contrasts
  psi_hat_ate <- psi_1_hat - psi_0_hat
  psi_hat_ate_direct <- mean(psi_ate_ic)

  # Point estimates + SEs relevant for binary outcomes
  Sig <- cov(cbind(psi_1_ic,psi_0_ic))
  psi_hat_rr <- psi_hat_or <- psi_hat_rr_se <- psi_hat_or_se <- NA

  if (y_bin) { # if outcome is binary
    psi_hat_rr <- psi_1_hat/psi_0_hat # risk ratio (only useful if binary)
    psi_hat_or <- (psi_1_hat/(1-psi_1_hat)) / (psi_0_hat/(1-psi_0_hat)) # odds ratio (only useful if binary)

    # standard errors (closed form via delta method)
    psi_hat_rr_se <- (Sig[1,1]/psi_0_hat^2 - 2*Sig[1,2]*psi_1_hat/(psi_0_hat^3) + Sig[2,2]*psi_1_hat^2/psi_0_hat^4)/n
    psi_hat_or_se <- psi_hat_or^2*(Sig[1,1]/(psi_1_hat^2*(1-psi_1_hat)^2) +
                                     Sig[2,2]/(psi_0_hat^2*(1-psi_0_hat)^2) -
                                     2*Sig[1,2]/(psi_1_hat*(1-psi_1_hat)*psi_0_hat*(1-psi_0_hat)))/n
  }

  return(list(ests = data.frame(psi_1_hat=mean(psi_1_ic),
                                psi_0_hat=mean(psi_0_ic),
                                psi_hat_ate=mean(psi_1_ic - psi_0_ic),
                                psi_hat_ate_direct=mean(psi_ate_ic),
                                psi_hat_rr=psi_hat_rr,
                                psi_hat_or=psi_hat_or
              ),
              vars = data.frame(psi_1_hat=var(psi_1_ic)/n,
                                psi_0_hat=var(psi_0_ic)/n,
                                psi_hat_ate=var(psi_1_ic - psi_0_ic)/n,
                                psi_hat_ate_direct=var(psi_ate_ic)/n,
                                psi_hat_rr=psi_hat_rr_se,
                                psi_hat_or=psi_hat_or_se
              ),
              ics = data.frame(psi_1_ic=psi_1_ic,
                               psi_0_ic=psi_0_ic,
                               psi_ate_ic=psi_ate_ic
              )
         )
  )
}
