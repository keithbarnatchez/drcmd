# drcmd

**D**oubly **R**obust **C**ausal Inference with **M**issing **D**ata

**Authors**: Keith Barnatchez and Griffin DesRoches

`drcmd` is an R package for implementing doubly-robust estimators of
counterfactual means in the presence of general missing data patterns.
`drcmd` leverages links between influence curves for counterfactual
means under no missingness, and the influence curve corresponding to the
missingness pattern in the user-supplied data. Detailed discussion of
the theoretical details behind the methods used in `drcmd` can be found
in Kennedy (2016), Tsiatis (2006), and van der Laan and Robins (2003).

Users can fit nuisance functions through Super Learner (a stacking
algorithm).

Please see the package vignette and documentation for further details.
Full documentation is available at the [package
website](https://kbarnatchez.com/drcmd/).

## Installation

``` r

devtools::install_github('keithbarnatchez/drcmd')
```

## Example

``` r

# Params for functions
eem_ind <- FALSE # TRUE = fit pseudo-outcome regression with empirical efficiency maximiztion
default_learners <- c('SL.glm','SL.gam') # default learners used for nuisance functions 
k <- 1 # number of cross-fitting folds
#-------------------------------------------------------------------------------
# Simulate simple missing outcome data structure
n <- 1e3
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X)) ; Y <- rnorm(n) + A + X
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) ; X <- as.data.frame(X)

# Make Y NA if R==0
Y[R==0] <- NA

drcmd_res <- drcmd::drcmd(Y,A,X, 
default_learners=default_learners,
eem_ind=eem_ind,k=k)
```

## Citation

Kennedy, E. H. (2016). *Semiparametric theory and empirical processes in
causal inference*. Statistical causal inferences and their applications
in public health research, 141-167.

Tsiatis, A. A. (2006). *Semiparametric theory and missing data* (Vol.
4). New York: Springer.

van der Laan, M. J., & Robins, J. M. (2003). *Unified methods for
censored longitudinal data and causality*. Springer New York.

## Contributing, reporting issues

- **Bugs/feature requests**: open an issue at
  <https://github.com/keithbarnatchez/drcmd/issues>  
- **Contributions**: PRs welcome — please run `devtools::check()` and
  `devtools::test()` before opening one  
- **Questions**: email <keithbarnatchez@gmail.com> or open a discussion
