# check_r_ind

Check if the missingness indicator R is defined correctly

## Usage

``` r
check_r_ind(data, Y, A, X, W, R)
```

## Arguments

- data:

  A data frame

- Y:

  A character string containing outcome variable name

- A:

  A character string containing treatment variable name

- X:

  A character vector containing covariate variable names

- W:

  A character vector containing weight variable names

- R:

  A character string containing randomization variable name

## Value

A logical value

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 200
R <- rbinom(n, 1, 0.7)
data <- data.frame(
  Y = ifelse(R == 1, rnorm(n), NA),
  A = ifelse(R == 1, rbinom(n, 1, 0.5), NA),
  X1 = rnorm(n),
  R = R
)
check_r_ind(data, Y = "Y", A = "A", X = "X1", W = character(0), R = "R")
} # }
```
