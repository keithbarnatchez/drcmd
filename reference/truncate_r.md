# Truncate complete case propensity scores

Truncate propensity scores to interval `[c, 1-c]`

## Usage

``` r
truncate_r(x, cutoff = 0.01)
```

## Arguments

- x:

  A vector of complete case propensity scores

## Value

A vector of complete propensity scores truncated to interval `[c, 1-c]`

## Examples

``` r
if (FALSE) { # \dontrun{
x <- c(0.001, 0.3, 0.5, 0.7, 0.999)
suppressWarnings(truncate_r(x, cutoff = 0.01))
} # }
```
