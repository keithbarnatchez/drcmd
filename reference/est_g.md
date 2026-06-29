# Estimate the propensity score

Function for obtaining propensity score estimates P(A=a\|X)

## Usage

``` r
est_g(idx, A, X, R, kappa_hat, g_learners, cv_folds = 5)
```

## Arguments

- idx:

  Indices to carry out estimation over

- A:

  A binary treatment variable (1=treated, 0=control)

- X:

  Dataframe containing baseline covariates

- R:

  Binary missingness indicator, where 0 indicates missing data

- kappa_hat:

  A numeric vector containing the fitted values of kappa

- g_learners:

  A character vector containing the names of the learners for estimation

- cv_folds:

  Number of cross-validation folds for SuperLearner

## Value

A list containing the estimate of `E[Y | A=a, X, W]`

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 500
X <- rnorm(n)
A <- rbinom(n, 1, plogis(X))
R <- rbinom(n, 1, plogis(X))
A[R == 0] <- 0
X <- data.frame(X)
kappa_hat <- rep(1, n)
g_learners <- c('SL.glm')
g_hat <- est_g(idx = 1:n, A = A, X = X, R = R,
               kappa_hat = kappa_hat, g_learners = g_learners)
} # }
```
