# Print drcmd object

S3 method for printing drcmd objects. Provides concise summary of
results, and prints values of optional arguments. Use summary() function
for more detailed summary of results.

## Usage

``` r
# S3 method for class 'drcmd'
print(x, ...)
```

## Arguments

- x:

  An object of class drcmd

- ...:

  Additional arguments (unused)

## Value

No return value. Called for printing a concise results summary

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
print(fit)
#> drcmd results
#> -------------------------------------------------
#> ATE estimate:  0.9474179 
#> -------------------------------------------------
#> Variables with missingness (U):  Y 
#> -------------------------------------------------
#> Variables without missingness (Z):  X A 
#> -------------------------------------------------
#> Validity of results requires causal assumptions to hold
#> As well as the assumption that U is independent of R given Z
#> Number of cross-fitting folds (k): 1 
```
