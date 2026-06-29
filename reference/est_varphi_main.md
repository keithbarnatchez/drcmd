# Perform pseudo-outcome regression

Function for obtaining estimate of `E[phi_a | Z]` through pseudo-outcome
regression. Calls inner function to perform estimation, depending on
whether user wishes to perform empirical efficiency maximization or not

## Usage

``` r
est_varphi_main(
  idx,
  R,
  Z,
  phi_1_hat,
  phi_0_hat,
  kappa_hat,
  eem_ind,
  po_learners,
  Y,
  cv_folds = 5,
  quiet = TRUE,
  att = FALSE,
  atc = FALSE,
  phi_att_hat = NULL,
  phi_atc_hat = NULL
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

  Vector of pseudo-outcome values under A=1

- phi_0_hat:

  Vector of pseudo-outcome values under A=0

- kappa_hat:

  A numeric vector containing the fitted values of kappa

- eem_ind:

  A logical value indicating whether to estimate via EEM

- po_learners:

  A character vector containing the names of the SuperLearner algorithms

- Y:

  Outcome variable. Can be continuous or binary

- cv_folds:

  Number of cross-validation folds for SuperLearner

- quiet:

  Logical indicating whether to suppress progress messages

- att:

  Logical indicating whether to compute pseudo-outcome regression for
  ATT

- atc:

  Logical indicating whether to compute pseudo-outcome regression for
  ATC

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
varphi <- est_varphi_main(1:n, R, Z, phi_1_hat, phi_0_hat,
                          kappa_hat, eem_ind = FALSE,
                          po_learners = "SL.glm", Y = Y)
} # }
```
