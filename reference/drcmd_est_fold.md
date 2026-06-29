# Calculate point estimates within a single cross-fitting fold

Outer function for obtaining point estimates for single fold

## Usage

``` r
drcmd_est_fold(
  splits,
  Y,
  A,
  X,
  Z,
  R,
  m_learners,
  g_learners,
  r_learners,
  po_learners,
  eem_ind,
  tml,
  Rprobs,
  cutoff,
  y_bin,
  cv_folds = 5,
  quiet = TRUE,
  att = FALSE,
  atc = FALSE
)
```

## Arguments

- splits:

  A list of train/test indices

- Y:

  Outcome variable. Can be continuous or binary

- A:

  A binary treatment variable (1=treated, 0=control)

- X:

  Dataframe containing baseline covariates

- Z:

  Dataframe containing all variables that are never subject to
  missingness

- R:

  Binary missingness indicator

- m_learners:

  A character vector containing learners to be used for the outcome
  regression. A vector of SuperLearner library names

- g_learners:

  A character vector containing learners to be used for the propensity
  score. A vector of SuperLearner library names

- r_learners:

  A character vector containing learners to be used for the missingness
  indicator regression. A vector of SuperLearner library names

- po_learners:

  A character vector containing learners to be used for the
  pseudo-outcome regression. A vector of SuperLearner library names

- eem_ind:

  A logical indicating whether to use empirical efficiency maximization

- Rprobs:

  A vector of probabilities for R

## Value

A list of results from estimation on current fold

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
splits <- list(train = 1:n, test = 1:n)
fold_res <- drcmd_est_fold(splits, Y, A, X, Z, R,
                           m_learners = "SL.glm", g_learners = "SL.glm",
                           r_learners = "SL.glm", po_learners = "SL.glm",
                           eem_ind = FALSE, tml = FALSE,
                           Rprobs = NA, cutoff = 0.025, y_bin = FALSE)
} # }
```
