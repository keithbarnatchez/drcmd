# dcrmd

**D**oubly **R**obust **C**ausal Inference with **M**issing **D**ata

**Author**: Keith Barnatchez

`drcmd` is an R package for implementing doubly-robust estimators of counterfactual means in the presence of general missing data patterns. `drcmd` leverages links between influence curves for counterfactual means under no missingngess, and the influence curve corresponding to the missingness pattern in the user-supplied data. Detailed discussion of -- can be found in Kennedy (2015) and Tisatis (2006).

Users can fit nuisance functions through either SuperLearner or the highly-adaptive LASSO. 

Please see the package vignette and documentation for further details.

------------------------------------------------------------------------
## Installation

```r
devtools::install_github('keithbarnatchez/drcmd')
```

------------------------------------------------------------------------
## Example

```r
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

------------------------------------------------------------------------
## Citation

------------------------------------------------------------------------
