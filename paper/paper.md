---
title: "`drcmd`: An R packagee for doubly-robust causal inference with missing data"
tags:
    - causal inference
    - missing data
    - R
authors: 
    - name: Keith Barnatchez
    - affiliation: 1
affiliations:
    - name: Harvard T. H. Chan School of Public Health
    - index: 1
date: 16 May 2025
bibliography: refs.bib
---

# Summary

The `drcmd` `R` package performs semi-parametric efficient estimation of causal effects of point exposures in settings where the data available to the researcher is subject to missingness. By implementing methods from semi-parametric theory and missing data (see, e.g. @robins1994estimation and @kennedy2022semiparametric), `drcmd` accomodates general patterns of missing data, while enabling users to estimate nuisance functions with flexible machine leaening methods. `drcmd` automatically determines the missingness patterns present in user-supplied data, and provides information on assumptions that must hold regarding the missingness mechanisms in order for point estimates and inferences to be valid. By accommodating *general* patterns of missingness, `drcmd` serves as a centralized library for researchers aiming to perform causal inference with missing data.

# Statement of Need
The use of doubly-robust methods for performing causal inference of point exposures on outcomes of interest has surged over the past decade, and numerous software packages have been developed for implementing these estimators. While these packages are well-suited for use on complete, missingness-free data, leading statistical software packages provide little to no support for missing data. The lack of a centralized package for performing doubly-robust causal inference has functioned as a severe impediment for researchers, as missing data is ubiquitous in real-world data, and the specific patterns of missingness can greatly vary across applications. `drcmd` addresses this shortcoming by providing a single R package for performing doubly-robust causal inference in the presence of general missing data patterns.

# Examples of Missing Data

To illustrate the breadth of missingness patterns `drcmd` can accommodate, we outline a non-exhaustive set of example problems that can be directly adressed by `drcmd`. Throughout, assume that the causal inference assumptions of consistency, positivity and unconfoundedness hold across all data structures provided below.

### Example 1: Two phase sampling 

The statistical theory `drcmd` leverages is most frequently invoked in two-phase sampling settings, where the researcher controls the proportion of missingness for variables that are typically difficult to measure. Suppose one observes 
$$
(Y_i,A_i, V_i, R_i U_i, R_i),  \ i = 1,\ldots,n
$$
where $Y$ is an outcome of interest, $A$ a binary treatment, $R$ a complete case indicator, and $X=(V,U)$ a vector of covariates sufficient to adjust for confounding. A subset of variables in $X$, denoted $U$, are expensive to measure, and only measured for the subset of units in which $R_i=1$. This subset is usually selected according to a sampling function dependent on the initially observed variables: $P(R=1| Y, A, V)$. 

Let $Y(a)$ denote a unit's counterfactual outcome under treatment level $A=a$. When interest lies in estimating general functions of counterfactual means---such as average treatment effects---`drcmd` can be readily applied to such a data structure. Further, the selection probabilities $P(R=1| Y, A, V)$ can be directly incorporated by the package.

### Example 2: Missing outcome and treatment

While the methods underpinning `drcmd` have most commonly leveraged in two-phase settings like the one above, their use is appliable in broader missing data problems. Suppose the observed data structure is
$$
(R_i Y_i, R_i A_i, X_i, R_i), \ i = 1,\ldots,n
$$
where $(Y, A) \perp R | X$. Unlike two-phase sampling problems, the missingness mechanism $P(R=1|X)$ must be estimated. `drcmd` readily handles such data structures.

### Example 3: Surrogate outcomes

In surrogate outcome settings, one observes
$$
(R_i Y_i, A_i, X_i, R_i, S_i), \ i = 1,\ldots,n
$$
where $S_i$ is an always-observed surrogate outcome that serves as a proxy for the partially missing outcome of interest $Y$. When $Y \perp R | X, A, S$, `drcmd` can produce valid treatment effect estimates.

# Availability 
`drcmd` is publically available for download on [GitHub](https://github.com/keithbarnatchez/drcmd). Further information on the use of the package can be found in the package vignette and user manual.

# Acknowledgements

This work was partially funded by National Institute of Health (NIH) grant T32AI007358.

# References

