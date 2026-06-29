# Cross-fit SE estimation via influence curve contributions

Estimates cross-fit SEs by taking contributions to influence curve of
each fold

## Usage

``` r
est_ses_crossfit(res, y_bin, att = FALSE, atc = FALSE)
```

## Arguments

- res:

  A list of per-fold results from drcmd_est_fold, each containing an
  `ics` data frame of influence curve contributions

- y_bin:

  Logical indicating whether the outcome is binary

- att:

  Logical indicating whether to compute SEs for the ATT

- atc:

  Logical indicating whether to compute SEs for the ATC

## Value

A list of cross-fit SEs

## Examples

``` r
if (FALSE) { # \dontrun{
# Typically called internally by drcmd_est() after cross-fitting.
# res is a list of per-fold results from drcmd_est_fold(), each
# containing an $ics data frame of influence curve contributions.
} # }
```
