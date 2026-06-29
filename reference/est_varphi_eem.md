# Perform pseudo-outcome regression with empirical efficiency maximization

Function for obtaining estimate of `E[phi_a | Z]` via empirical
efficiency maximization

## Usage

``` r
est_varphi_eem(
  idx,
  R,
  Z,
  phi_1_hat,
  phi_0_hat,
  kappa_hat,
  po_learners,
  Y,
  cv_folds = 5
)
```

## Arguments

- idx:

  Indices to carry out estimation over

- R:

  Binary missingness indicator, where 0 indicates missing data

- Z:

  Dataframe containing all non-missing variables

- phi_1_hat:

  Vector of predicted phi under A=1

- phi_0_hat:

  Vector of predicted phi under A=0

- kappa_hat:

  A numeric vector containing the fitted values of kappa

- po_learners:

  A character vector containing the names of the superlearner algorithms

- Y:

  Outcome variable. Can be continuous or binary

## Value

A list containing the estimate of `E[phi_a | Z]` for a=0 and a=1

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 500
Z <- data.frame(Z1 = rnorm(n))
R <- rbinom(n, 1, 0.7)
phi_1_hat <- rnorm(n)
phi_0_hat <- rnorm(n)
kappa_hat <- rep(0.7, n)
Y <- rnorm(n)
varphi <- est_varphi_eem(1:n, R, Z, phi_1_hat, phi_0_hat,
                         kappa_hat, po_learners = "SL.glm", Y = Y)
} # }
```
