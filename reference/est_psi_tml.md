# Estimate causal effects through targeted maximum likelihood

Function for obtaining estimates of causal estimands through targeted
maximum likelihood. To be called within drcmd_est_fold(). Runs a
separate TML fluctuation for each of the ATE, `E[Y(1)]` and `E[Y(0)]`,
so each estimand is debiased directly. ATT/ATC are not supported here –
use the one-step option (tml=FALSE) for those.

## Usage

``` r
est_psi_tml(
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
  varphi_hat,
  att = FALSE,
  atc = FALSE
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

- att:

  Logical. Not supported under TML; must be FALSE

- atc:

  Logical. Not supported under TML; must be FALSE

## Value

A list with point estimates, variances, and influence-curve
contributions for `E[Y(1)]`, `E[Y(0)]`, and ATE

## Examples

``` r
if (FALSE) { # \dontrun{
# Typically called internally by drcmd_est_fold() when tml = TRUE.
# Requires nuisance estimates from get_nuisance_ests() and
# pseudo-outcome estimates from est_varphi_main().
} # }
```
