# drcmd

**D**oubly **R**obust **C**ausal Inference with **M**issing **D**ata

[![R-CMD-check](https://github.com/keithbarnatchez/drcmd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/keithbarnatchez/drcmd/actions/workflows/R-CMD-check.yaml)

**Authors**: Keith Barnatchez and Griffin DesRoches

`drcmd` is an R package for implementing doubly-robust estimators of counterfactual means in the presence of general missing data patterns. `drcmd` leverages links between influence curves for counterfactual means under no missingness, and the influence curve corresponding to the missingness pattern in the user-supplied data. Detailed discussion of the theoretical details behind the methods used in `drcmd` can be found in Kennedy (2016), Tsiatis (2006), and van der Laan and Robins (2003).

Users can fit nuisance functions through Super Learner (a stacking algorithm). 

Please see the package vignette and documentation for further details. Full documentation is available at the [package website](https://kbarnatchez.com/drcmd/).

## Installation

Install the development version from GitHub with `remotes`:

```r
install.packages("remotes")
remotes::install_github("keithbarnatchez/drcmd")
```

Core dependencies are installed automatically. Some optional Super Learner
libraries require suggested packages such as `gam`, `hal9001`, or `nnls`; install
the package corresponding to any optional learner you include.

## Example

```r
set.seed(1)

# Simulate a simple missing-outcome setting
n <- 1000
X <- rnorm(n)
A <- rbinom(n, 1, plogis(X / 2))
Y <- rnorm(n) + A + X
Ystar <- Y + rnorm(n) / 2
R <- rbinom(n, 1, plogis(X / 2))

# The outcome is observed only for complete cases
Y[R == 0] <- NA

fit <- drcmd::drcmd(
  Y = Y,
  A = A,
  X = data.frame(X = X),
  W = data.frame(Ystar = Ystar),
  default_learners = c("SL.glm", "SL.gam"),
  eem_ind = FALSE,
  k = 1 # Set greater than 1 to enable cross-fitting
)

summary(fit)
```

## Citation

If you use `drcmd` in your work, please cite the package:

Barnatchez, K. and DesRoches, G. (2026). *drcmd: Doubly-Robust Causal Inference with Missing Data*. R package version 0.1.0. https://github.com/keithbarnatchez/drcmd

A BibTeX entry is available from R:

```r
citation('drcmd')
```

## References

Kennedy, E. H. (2016). *Semiparametric theory and empirical processes in causal inference*. Statistical causal inferences and their applications in public health research, 141-167.

Tsiatis, A. A. (2006). *Semiparametric theory and missing data* (Vol. 4). New York: Springer.

van der Laan, M. J., & Robins, J. M. (2003). *Unified methods for censored longitudinal data and causality*. Springer New York.

## Contributing, reporting issues

- **Bugs and feature requests**: open an issue at https://github.com/keithbarnatchez/drcmd/issues
- **Contributions**: pull requests are welcome; please run `devtools::check()` and `devtools::test()` before opening one
- **Questions and support**: email keithbarnatchez@gmail.com
