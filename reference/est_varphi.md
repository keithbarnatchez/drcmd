# Perform pseudo-outcome regression with conventional loss function

Outer function for obtaining estimate of `E[phi | Z]`

## Usage

``` r
est_varphi(idx, R, Z, phi_1_hat, phi_0_hat, po_learners, Y, cv_folds = 5)
```

## Arguments

- idx:

  Indices to carry out estimation over

- R:

  Binary missingness indicator, where 0 indicates missing data

- Z:

  Dataframe containing all non-missing variables

- phi_1_hat:

  A numeric vector containing the fitted values of phi under A=1

- phi_0_hat:

  A numeric vector containing the fitted values of phi under A=0

- po_learners:

  A character vector containing the names of the superlearner algorithms

## Value

A list containing the estimate of `E[phi | Z]` for a=0 and a=1

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 500
Z <- data.frame(Z1 = rnorm(n))
R <- rbinom(n, 1, 0.7)
phi_1_hat <- rnorm(n)
phi_0_hat <- rnorm(n)
Y <- rnorm(n)
varphi <- est_varphi(1:n, R, Z, phi_1_hat, phi_0_hat,
                     po_learners = "SL.glm", Y = Y)
} # }
```
