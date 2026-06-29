# Estimating outcome regression

Function for obtaining estimate of `E[Y | A=a, X]`. Estimation is
carried out with SuperLearer, using learners specified in m_learners

## Usage

``` r
est_m_a(idx, Y, A, X, R, kappa_hat, m_learners, cv_folds = 5)
```

## Arguments

- idx:

  Indices to carry out estimation over

- Y:

  Outcome variable. Can be continuous or binary

- A:

  A binary treatment variable (1=treated, 0=control)

- X:

  Dataframe containing baseline covariates

- R:

  Binary missingness indicator, where 0 indicates missing data

- kappa_hat:

  A numeric vector containing the fitted values of kappa

- m_learners:

  A character vector containing the names of the superlearner algorithms

## Value

A list containing the estimate of `E[Y | A=a, X]`

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 500
X <- rnorm(n)
A <- rbinom(n, 1, plogis(X))
R <- rbinom(n, 1, plogis(X))
Y <- A + X + rnorm(n)
Y[R == 0] <- 0; A[R == 0] <- 0
X <- data.frame(X)
kappa_hat <- rep(1, n)
m_learners <- c('SL.glm')
m_hat <- est_m_a(idx = 1:n, Y = Y, A = A, X = X, R = R,
                 kappa_hat = kappa_hat, m_learners = m_learners)
} # }
```
