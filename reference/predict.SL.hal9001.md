# Prediction wrapper for the highly-adaptive lasso

Prediction wrapper for the highly-adaptive lasso (HAL) algorithm
implemented in the hal9001 package

## Usage

``` r
# S3 method for class 'SL.hal9001'
predict(object, newdata, ...)
```

## Arguments

- object:

  A fitted SL.hal9001 object

- newdata:

  Data frame of covariates for prediction

- ...:

  Additional arguments (unused)

## Value

Predicted values on new data

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 200
X <- data.frame(X1 = rnorm(n))
Y <- X$X1 + rnorm(n)
fit <- SL.hal9001(Y, X, newX = X, family = gaussian(),
                  obsWeights = rep(1, n))
preds <- predict(fit$fit, newdata = X)
} # }
```
