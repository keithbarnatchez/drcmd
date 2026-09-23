# Plot results from drcmd object

S3 method for plotting results from drcmd object. Plots are available
for the following: (1) pseudo-outcome regression fit, (2) influence
curve distribution, (3) treatment propensity score distribution, and (4)
complete-case propensity score distribution. When type='All' (the
default), user can view all four plots in succession interactively

## Usage

``` r
# S3 method for class 'drcmd'
plot(x, type = "All", ...)
```

## Arguments

- x:

  An object of class drcmd

- type:

  Character denoting type of plot to generate. Must be one of 'All',
  'PO', 'IC', 'g_hat', 'r_hat'

- ...:

  Additional arguments (unused)

## Value

No return value. Called for plotting results from drcmd object

## Examples

``` r
set.seed(1)
n <- 200
X <- rnorm(n)
A <- rbinom(n, 1, plogis(X))
Y <- rnorm(n) + A + X
R_ind <- rbinom(n, 1, plogis(X))
Y[R_ind == 0] <- NA
covariates <- data.frame(X = X)

fit <- drcmd(Y, A, covariates, default_learners = "SL.glm", k = 1)
#> Warning: Augmentation regression selected an all-zero ensemble; its predictions are zero. Consider revising the learner library.
plot(fit, type = "PO")
```
