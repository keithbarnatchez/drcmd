# Doubly-robust causal inference with missing data

Doubly-robust estimation of counterfactual means and average treatment
effects for binary point treatments when outcomes, treatments, or
covariates are incompletely observed, with risk ratios and odds ratios
available for binary outcomes. Identification relies on consistency,
conditional treatment exchangeability given `X`, treatment positivity,
and independence of the complete-case indicator and partially observed
variables conditional on fully observed variables, with positive
complete-case probabilities. Missing-data methods follow the
semiparametric framework described by Tsiatis (2006). Nuisance functions
are estimated using 'SuperLearner' with user-specified libraries and
optional cross-fitting. Known complete-case probabilities can be
supplied through `Rprobs` for two-phase sampling designs, and fully
observed auxiliary variables can be included in the missingness and
augmentation regressions through `W`.

## Usage

``` r
drcmd(
  Y,
  A,
  X,
  W = NA,
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

  (optional) Dataframe containing fully observed auxiliary variables,
  such as proxy measurements of partially observed variables. These
  variables enter the missingness and augmentation regressions, but not
  the outcome or treatment regressions conditional on `X`; covariates
  needed for treatment exchangeability should therefore be included in
  `X`.

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

## Details

Treatment and complete-case regressions stop with an error if the fitted
ensemble has only zero coefficients or predicts only zeros or only ones.
Such fits are rejected before probability truncation. A zero-valued
augmentation fit is permitted, but emits a `drcmd_zero_augmentation`
warning so that the learner specification can be reviewed.

TML omits the missingness update when all observations are complete and
complete-case probabilities equal one. An update with a zero clever
covariate on all informative evaluation observations is also left
unchanged. Other targeting fits must converge to an identifiable finite
coefficient; otherwise estimation stops. Outcome targeting uses a
quasibinomial logit regression, which retains the binomial estimating
equation for fractional outcomes and inverse-probability weights without
interpreting them as binomial counts.

## References

Tsiatis, A. A. (2006). Semiparametric Theory and Missing Data. Springer.
[doi:10.1007/0-387-37345-4](https://doi.org/10.1007/0-387-37345-4) .

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
#> Warning: Augmentation regression selected an all-zero ensemble; its predictions are zero. Consider revising the learner library.
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
