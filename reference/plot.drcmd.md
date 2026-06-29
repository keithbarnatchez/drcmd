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
if (FALSE) { # \dontrun{
n <- 2500
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <-  rbinom(n,1,plogis(X-A)) # rnorm(n) + A + X + X^2 + A*X + sin(X) # note: true ATE is 1
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) # error-prone outcome measurements

# Make A NA if R==0
A[R==0] <- NA

# Obtain ATE estimates, fitting all nuisance models with ensemble of splines +
# GAMs (save for the pseudo-outcome regression, which is done with XGboost)
drcmd_res <- drcmd(Y,A,covariates,
                   default_learners= c('SL.gam','SL.glm'),
                   po_learners = 'SL.gam')
plot(drcmd_res,type='PO')
} # }
```
