# Contruct full data EIF estimate

Function for constructing estimate of the full data EIF, given
propensity score and outcome model fits

## Usage

``` r
get_phi_hat(Y, A, X, R, Z, g_hat, m_a_hat, kappa_hat)
```

## Arguments

- Y:

  Outcome variable. Can be continuous or binary

- A:

  A binary treatment variable (1=treated, 0=control)

- X:

  Dataframe containing baseline covariates

- R:

  Binary missingness indicator

- Z:

  Dataframe containing all variables that are never subject to
  missingness

- g_hat:

  Propensity score predictions

- m_a_hat:

  List of outcome predictions under A=0/1

- kappa_hat:

  Missingness probabilities

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 500
X <- data.frame(X1 = rnorm(n))
A <- rbinom(n, 1, plogis(X$X1))
Y <- A + X$X1 + rnorm(n)
R <- rep(1, n)
Z <- X
g_hat <- plogis(X$X1)
m_a_hat <- list(m_1_hat = 1 + X$X1, m_0_hat = X$X1)
kappa_hat <- rep(1, n)
phi <- get_phi_hat(Y, A, X, R, Z, g_hat, m_a_hat, kappa_hat)
} # }
```
