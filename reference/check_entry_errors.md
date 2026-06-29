# Check arguments to drcmd for entry errors

Checks if the arguments to main function drcmd are correctly specified.
Returns TRUE if all checks are passed, throws an error with message
otherwise.

## Usage

``` r
check_entry_errors(Y, A, X, W, R, eem_ind, Rprobs, k)
```

## Arguments

- Y:

  A vector or data frame containing outcome values

- A:

  A vector or data frame containing treatment variable values

- X:

  A data frame containing covariate values

- W:

  A data frame containing proxy variable values

- R:

  A vector containing missingness indicator variable

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 200
X <- data.frame(X1 = rnorm(n))
W <- data.frame(W1 = rnorm(n))
A <- rbinom(n, 1, 0.5)
Y <- rnorm(n)
R <- rbinom(n, 1, 0.7)
check_entry_errors(Y, A, X, W, R, eem_ind = FALSE, Rprobs = NA, k = 1)
} # }
```
