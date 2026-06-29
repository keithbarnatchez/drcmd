# Doubly-robust causal inference with missing data

Doubly-robust estimation of counterfactual means in the presence of
missing data. The drcmd() function estimates counterfactual means for
binary point treatments and reports average treatment effects, as well
as causal risk ratios and odds ratios for binary outcomes. General
missingness patterns in the data are allowed and automatically
determined by the function – the only requirement is that any
missingness occurs at random conditional on variables that are always
available. For scenarios where non-missingness probabilities are known,
as is common in two-phase sampling designs, users can provide the
non-missingness probabilities through the Rprobs argument. Users can fit
nuisance functions through either highly-adaptive LASSO (HAL) or
SuperLearner, the latter of which the user must specify libraries.

## Usage

``` r
drcmd(
  Y,
  A,
  X,
  W = NA,
  R = NA,
  default_learners = NULL,
  m_learners = NULL,
  g_learners = NULL,
  r_learners = NULL,
  po_learners = NULL,
  eem_ind = FALSE,
  tml = FALSE,
  Rprobs = NA,
  k = 1,
  cutoff = 0.025,
  quiet = TRUE,
  cv_folds = 5,
  parallel = FALSE,
  att = FALSE,
  atc = FALSE
)
```

## Arguments

- Y:

  Outcome variable. Can be continuous or binary

- A:

  A binary treatment variable (1=treated, 0=control)

- X:

  Dataframe containing baseline covariates

- W:

  (optional) Dataframe containing variables solely predictive of
  missingness, but not a cause of the outcome or exposure.

- R:

  (optional) A character string specifying the missingness indicator,
  where 0 indicates missing data. If not specified, the function will
  create the missingness indicator by identifying the missingness
  pattern in the data

- default_learners:

  A character vector containing SuperLearner libraries to use for
  estimating all nuisance functions. User can alternatively specify
  libraries for each nuisance function for added flexibility

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

- tml:

  A logical indicating whether to obtain estimates through targeted
  maximum likelihood estimation (TML). If TRUE, estimates obtained by
  TML. If FALSE (default setting), estimates are obtained through a
  one-step estimator

- Rprobs:

  A vector of probabilities for the missingness indicator. Only suitable
  for study designs where researcher controls mechanism by which
  variables are missing (e.g. two-phase sample designs). Defaults to NA,
  in which case missingness probabilities are estimated.

- k:

  A numeric indicating the number of folds for cross-fitting

- cutoff:

  Cutoff for treatment and complete case propensity scores. Estimates
  outside of `[cutoff, 1-cutoff]` are set to cutoff or 1-cutoff,
  respectively

- quiet:

  Logical indicating whether to suppress progress messages. Default is
  TRUE. Set to FALSE to see which nuisance function is being estimated.

- cv_folds:

  Number of cross-validation folds used internally by SuperLearner for
  model selection. Default is 5. Lower values speed up estimation.

- parallel:

  Logical indicating whether to run cross-fitting folds in parallel
  using
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html). Only
  used when k \> 1. Note: not supported on Windows.

- att:

  Logical indicating whether to additionally estimate the average
  treatment effect on the treated (ATT). Default is FALSE.

- atc:

  Logical indicating whether to additionally estimate the average
  treatment effect on the controls (ATC). Default is FALSE.

## Value

An S3 object of class `"drcmd"` containing estimation results,
information on the missing data structure, and parameters used in the
estimation

- `params`:

  A list containing parameters used in the estimation

- `Z`:

  A character vector containing always-available variables

- `R`:

  A character vector containing the complete case indicator values

- `U`:

  A character vector containing partially-missing variables

- `results`:

  A list of dataframes storing (i) point estimates, (ii) standard
  errors, and (iii) nuisance function estimates

## Examples

``` r
set.seed(1)
n <- 200
X <- rnorm(n)
A <- rbinom(n, 1, plogis(X))
Y <- rnorm(n) + A + X
R_ind <- rbinom(n, 1, plogis(X))
Y[R_ind == 0] <- NA
covariates <- data.frame(X = X)

fit <- drcmd(Y, A, covariates,
             default_learners = "SL.glm",
             k = 1)
#> Loading required package: nnls
fit
#> drcmd results
#> -------------------------------------------------
#> ATE estimate:  0.9474179 
#> -------------------------------------------------
#> Variables with missingness (U):  Y 
#> -------------------------------------------------
#> Variables without missingness (Z):  X A 
#> -------------------------------------------------
#> Validity of results requires causal assumptions to hold
#> As well as the assumption that U is independent of R given Z
#> Number of cross-fitting folds (k): 1 
```
