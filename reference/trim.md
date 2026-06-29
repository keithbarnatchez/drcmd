# Trim vector (for numerical stability)

Trims values of a vector to avoid numerical instability

## Usage

``` r
trim(x, val = .Machine$double.neg.eps)
```

## Arguments

- x:

  A vector of values

- val:

  A small value to add to 0 and subtract from 1

## Value

A vector with values trimmed to avoid numerical instability
