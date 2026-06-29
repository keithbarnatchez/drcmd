# Obtain nuisance function estimates

Function for obtaining estimates of all relevant nuisance functions

## Usage

``` r
get_nuisance_ests(
  idx,
  Y,
  A,
  X,
  Z,
  R,
  m_learners,
  g_learners,
  r_learners,
  Rprobs,
  cutoff,
  cv_folds = 5,
  quiet = TRUE
)
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

- Z:

  Dataframe containing all non-missing variables

- R:

  Binary missingness indicator, where 0 indicates missing data

- m_learners:

  A character vector containing learners to be used for the outcome
  regression. A vector of SuperLearner library names

- g_learners:

  A character vector containing learners to be used for the propensity
  score. A vector of SuperLearner library names

- r_learners:

  A character vector containing learners to be used for the missingness
  indicator regression. A vector of SuperLearner library names

- Rprobs:

  A vector of probabilities for R

- cutoff:

  Cutoff for propensity score truncation

- cv_folds:

  Number of cross-validation folds for SuperLearner

- quiet:

  Logical indicating whether to suppress progress messages

## Value

A list of nuisance estimates

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 500
X <- data.frame(X1 = rnorm(n))
A <- rbinom(n, 1, plogis(X$X1))
Y <- A + X$X1 + rnorm(n)
R <- rbinom(n, 1, plogis(X$X1))
Y[R == 0] <- 0; A[R == 0] <- 0
Z <- X
nuis <- get_nuisance_ests(1:n, Y, A, X, Z, R,
                          m_learners = "SL.glm",
                          g_learners = "SL.glm",
                          r_learners = "SL.glm",
                          Rprobs = NA, cutoff = 0.025)
} # }
```
