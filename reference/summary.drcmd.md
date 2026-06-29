# Summarize results from drcmd

S3 method for summarizng drcmd results. Provides detailed summary of
estimation output, while also providing information on missingness
mechanism. Can optionally print out values of user-supplied arguments by
setting detail=TRUE

## Usage

``` r
# S3 method for class 'drcmd'
summary(object, detail = FALSE, ...)
```

## Arguments

- object:

  An object of class drcmd

- detail:

  Logical. If TRUE, print out values of user-supplied arguments

- ...:

  Additional arguments (unused)

## Value

No return value. Called for printing a detailed results summary

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
summary(fit)
#> ======================================================================
#>                         Summary of drcmd results                      
#> ======================================================================
#> Estimand            Estimate           SE                     95% CI
#> ----------------------------------------------------------------------
#> ATE:                   0.947        0.222             [0.513, 1.382]
#> E[Y(1)]:               0.853        0.175             [0.511, 1.196]
#> E[Y(0)]:              -0.094        0.187            [-0.460, 0.272]
#> ----------------------------------------------------------------------
#> Variables with missingness (U): Y
#> Variables without missingness (Z): X, A
#> ----------------------------------------------------------------------
#> Validity of results requires causal assumptions to hold,
#> as well as the assumption that U is independent of R given Z
```
