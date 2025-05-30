% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drcmd.R
\name{drcmd_est_fold}
\alias{drcmd_est_fold}
\title{Calculate point estimates within a single cross-fitting fold}
\usage{
drcmd_est_fold(
  splits,
  Y,
  A,
  X,
  Z,
  R,
  m_learners,
  g_learners,
  r_learners,
  po_learners,
  eem_ind,
  tml,
  Rprobs,
  cutoff,
  y_bin
)
}
\arguments{
\item{splits}{A list of train/test indices}

\item{Y}{Outcome variable. Can be continuous or binary}

\item{A}{A binary treatment variable (1=treated, 0=control)}

\item{X}{Dataframe containing baseline covariates}

\item{Z}{Dataframe containing all variables that are never subject to missingness}

\item{R}{Binary missingness indicator}

\item{m_learners}{A character vector containing learners to be used for the
outcome regression. User can specify 'hal' or a vector of SuperLearner libraries}

\item{g_learners}{A character vector containing learners to be used for the
propensity score. User can specify 'hal' or a vector of SuperLearner libraries}

\item{r_learners}{A character vector containing learners to be used for the
missingness indicator regression. User can specify 'hal' or a vector of
SuperLearner libraries}

\item{po_learners}{A character vector containing learners to be used for the
pseudo-outcome regression. User can specify 'hal' or a vector of SuperLearner libraries}

\item{eem_ind}{A logical indicating whether to use empirical efficiency maximization}

\item{Rprobs}{A vector of probabilities for R}

\item{c}{}
}
\value{
A list of results from estimation on current fold
}
\description{
Outer function for obtaining point estimates for single fold
}
