# SuperLearner wrapper for the highly-adaptive lasso

Wrapper for the highly-adaptive lasso (HAL) algorithm implemented in the
hal9001 package

## Usage

``` r
SL.hal9001(Y, X, newX, family, obsWeights, ...)
```

## Arguments

- Y:

  Outcome variable

- X:

  Matrix of covariates

- newX:

  Matrix of covariates for prediction

- family:

  Family of the outcome variable

- obsWeights:

  Observation weights

- ...:

  Additional arguments to pass to the HAL function

## Value

A list containing the prediction and the fitted model

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 200
X <- data.frame(X1 = rnorm(n))
Y <- X$X1 + rnorm(n)
fit <- SL.hal9001(Y, X, newX = X, family = gaussian(),
                  obsWeights = rep(1, n))
head(fit$pred)
} # }
```
