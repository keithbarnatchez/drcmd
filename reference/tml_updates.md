# Update nuisance functions with TML

Function for updating the nuisance functions m_a and kappa through TML

## Usage

``` r
tml_updates(
  idx,
  Y,
  A,
  X,
  R,
  Z,
  m_1_hat,
  m_0_hat,
  g_hat,
  kappa_hat,
  phi_1_hat,
  phi_0_hat,
  varphi_hat
)
```

## Arguments

- idx:

  Indices of the training set

- Y:

  Outcome vector

- A:

  Treatment vector

- X:

  Dataframe of covariates

- R:

  Complete case indicator vector

- Z:

  Dataframe containing all variables that are never subject to
  missingness

- m_1_hat:

  A numeric vector containing the fitted outcome predictions under A=1

- m_0_hat:

  A numeric vector containing the fitted outcome predictions under A=0

- g_hat:

  A numeric vector containing the fitted treatment propensity scores

- kappa_hat:

  A numeric vector containing the fitted values of kappa

- phi_1_hat:

  Estimates of EIC for `E[Y(1)]`

- phi_0_hat:

  Estimates of `E[Y(0)]`

- varphi_hat:

  A list containing estimate of `E[phi_a | Z]`

## Value

A list containing updated estimates of the nuisance parameters

## Examples

``` r
if (FALSE) { # \dontrun{
# Typically called internally by est_psi_tml() to perform the TML
# targeting step on outcome and missingness nuisance functions.
} # }
```
