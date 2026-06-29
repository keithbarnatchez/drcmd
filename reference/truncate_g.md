# Truncate treatment propensity scores

Truncate propensity scores to interval `[c, 1-c]`

## Usage

``` r
truncate_g(x, cutoff = 0.025)
```

## Arguments

- x:

  A vector of treatment propensity scores

## Value

A vector of treatment propensity scores truncated to interval `[c, 1-c]`

## Examples

``` r
if (FALSE) { # \dontrun{
x <- c(0.001, 0.3, 0.5, 0.7, 0.999)
suppressWarnings(truncate_g(x, cutoff = 0.025))
} # }
```
