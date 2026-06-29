# Obtain estimates for current cross fitting fold

Function for obtaining counterfactual mean estimates for the current
fold of the cross-fitting procedure

## Usage

``` r
est_psi(
  idx,
  R,
  Z,
  kappa_hat,
  phi_hat,
  varphi_hat,
  y_bin,
  att = FALSE,
  atc = FALSE,
  plugin_att = NULL,
  plugin_atc = NULL,
  phi_att_hat = NULL,
  phi_atc_hat = NULL
)
```

## Arguments

- idx:

  A data frame

- R:

  A character string containing randomization variable name

- Z:

  A character vector containing the names of the variables in Z

- kappa_hat:

  A numeric vector containing the fitted values of kappa

- phi_hat:

  A list containing the estimate of `E[phi | Z]`

- varphi_hat:

  A list containing the estimate of `E[phi | Z]`

## Value

A list of point estimates and standard errors for counterfactual means
and various counterfactual contrasts

## Examples

``` r
if (FALSE) { # \dontrun{
# Typically called internally by drcmd_est_fold() after obtaining
# nuisance estimates, phi_hat, and varphi_hat.
} # }
```
