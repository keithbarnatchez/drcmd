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

## Please see the package vignette and documentation for further details. Full documentation is available at the [package website](https://keithbarnatchez.github.io/drcmd/).

## Installation

``` r

devtools::install_github('keithbarnatchez/drcmd')
```

| \## Example |
|:---|
| \## Citation |
| Kennedy, E. H. (2016). *Semiparametric theory and empirical processes in causal inference*. Statistical causal inferences and their applications in public health research, 141-167. |
| Tsiatis, A. A. (2006). *Semiparametric theory and missing data* (Vol. 4). New York: Springer. |
| van der Laan, M. J., & Robins, J. M. (2003). *Unified methods for censored longitudinal data and causality*. Springer New York. |

## Contributing, reporting issues

- **Bugs/feature requests**: open an issue at
  <https://github.com/keithbarnatchez/drcmd/issues>  
- **Contributions**: PRs welcome — please run `devtools::check()` and
  `devtools::test()` before opening one  
- **Questions**: email <keithbarnatchez@gmail.com> or open a discussion
