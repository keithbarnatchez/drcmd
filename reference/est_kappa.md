# Estimate complete case propensity scores

Function for obtaining complete case propensity score estimates

## Usage

``` r
est_kappa(idx, Z, R, r_learners, cv_folds = 5)
```

## Arguments

- idx:

  Indices to carry out estimation over

- Z:

  Dataframe containing all non-missing variables

- R:

  Binary missingness indicator, where 0 indicates missing data

- r_learners:

  A character vector specifying learners to be used for estimation

## Value

A numeric vector containing the estimate of `E[R | Z]`

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 500
Z <- rnorm(n)
R <- rbinom(n, 1, plogis(Z))
Z <- data.frame(Z)
r_learners <- c('SL.glm')
kappa_hat <- est_kappa(idx = 1:n, Z = Z, R = R,
                       r_learners = r_learners)
} # }
```
