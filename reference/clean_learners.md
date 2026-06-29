# Clean SuperLearner libraries

Internal function for setting SuperLearner libraries before they are
passed into estimation procedures. Default libraries from the
sl_learners argument are used to fill in missing libraries for any
nuisance function with libraries left unspecified

## Usage

``` r
clean_learners(
  default_learners,
  m_learners,
  g_learners,
  r_learners,
  po_learners
)
```

## Arguments

- default_learners:

  Either null, or a character vector containing SuperLearner libraries
  to use for estimating all nuisance functions. User can alternatively
  specify libraries for each nuisance function for added flexibility

- m_learners:

  Either null, or a character vector containing SuperLearner libraries
  to be used for the outcome regression

- g_learners:

  Either null, or a character vector containing SuperLearner libraries
  to be used for the propensity scores

- r_learners:

  Either null, or a character vector containing SuperLearner libraries
  to be used for the missingness indicator regression

- po_learners:

  Either null, or a character vector containing SuperLearner libraries
  to be used for the pseudo outcome regression

## Value

A list of
