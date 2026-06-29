# find_missing_pattern

Find the missing pattern in the data

## Usage

``` r
find_missing_pattern(Y, A, X, W)
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

## Value

A character string containing the missing pattern

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 200
X <- data.frame(X1 = rnorm(n))
A <- rbinom(n, 1, 0.5)
Y <- rnorm(n)
R <- rbinom(n, 1, 0.7)
Y[R == 0] <- NA
W <- X[, 0]
result <- find_missing_pattern(Y, A, X, W)
result$Z  # variables without missingness
result$U  # variables with missingness
result$R  # complete case indicator
} # }
```
