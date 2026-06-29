---
title: "`drcmd`: An R package for doubly-robust causal inference with missing data"
tags:
  - causal inference
  - missing data
  - R
authors:
  - name: Keith Barnatchez
    orcid: 0000-0002-1451-9335
    affiliation: 1
  - name: Griffin DesRoches
    affiliation: 2
affiliations:
  - name: Johns Hopkins Bloomberg School of Public Health
    index: 1
  - name: University of Massachusetts Amherst
    index: 2
date: 16 June 2026
bibliography: references.bib
---

# Summary

The `drcmd` `R` package performs semi-parametric efficient estimation of causal effects of point exposures in settings where the data available to the researcher is subject to missingness. By implementing methods from semi-parametric theory and missing data (see, e.g. @robins1994estimation and @kennedy2022semiparametric), `drcmd` accommodates general patterns of missing data, while enabling users to estimate nuisance functions with flexible machine learning methods. `drcmd` automatically determines the missingness patterns present in user-supplied data, and provides information on assumptions that must hold regarding the missingness mechanisms in order for point estimates and inferences to be valid. By accommodating *general* patterns of missingness, `drcmd` serves as a centralized library for researchers aiming to perform causal inference with missing data.

# Statement of Need
The use of doubly-robust methods for performing causal inference of point exposures on outcomes of interest has surged over the past decade, and numerous software packages have been developed for implementing these estimators. While these packages are well-suited for use on complete, missingness-free data, leading statistical software packages provide incomplete support for missing data. The lack of a centralized package for performing doubly-robust causal inference has been a significant drawback for researchers, as missing data is ubiquitous in real-world data, and the specific patterns of missingness can greatly vary across applications. 

# State of the Field

A number of `R` packages implement doubly-robust estimators of causal effects, including `AIPW` [@zhong2021aipw], `tmle` [@gruber2012tmle], and `drtmle` [@benkeser2017doubly]. These packages are designed primarily for complete-data settings, and their support for missing data is limited to specific, isolated patterns. `AIPW` and `tmle` accommodate missingness only in the outcome, while `drtmle` additionally allows missingness in the treatment. None of these packages address missing covariates or arbitrary combinations of missing variables. The `twoStageDesignTMLE` package handles missing covariates, focusing on data obtained from two-phase sampling designs [@rose2011targeted], but similarly does not accomodate missingness in other variables. 

`drcmd` unifies these capabilities under a single interface. Its central observation is that the semi-parametric machinery developed for two-phase sampling---where the probability of observing expensive-to-measure variables is fixed by design---applies far more broadly to missing-data problems in which the missingness mechanism must instead be estimated [@robins1994estimation; @kennedy2022semiparametric]. Two-phase sampling emerges as a special case in which the non-missingness probabilities are known rather than estimated. By automatically detecting the missingness pattern present in user-supplied data, `drcmd` accommodates missingness in any combination of outcomes, treatments, and covariates within a single framework.

# Software Design

`drcmd` is organized around a single user-facing function, `drcmd()`, which takes an outcome `Y`, a binary treatment `A`, and a data frame of covariates `X`, with missing values supplied as `NA`. Variables predictive of missingness alone may optionally be passed through `W`. The package automatically detects the missingness pattern, partitioning the data into always-observed variables, a complete-case indicator, and partially-missing variables, and reports the conditional independence assumption required for valid inference.

Estimation proceeds by constructing the observed-data efficient influence function from its full-data counterpart [@tsiatis2006semiparametric]. `drcmd` fits four nuisance functions---an outcome regression, a treatment propensity score, a missingness mechanism, and a pseudo-outcome regression that projects the full-data influence function onto the always-observed variables---using SuperLearner with user-specified libraries [@vanderlaan2007super] and cross-fitting. Point estimates are obtained through either a one-step correction (the default) or targeted maximum likelihood estimation. When non-missingness probabilities are known by design, as in two-phase sampling, they may be supplied through `Rprobs` rather than estimated.

`drcmd()` returns an S3 object with `print()`, `summary()`, and `plot()` methods. The function reports counterfactual means, average treatment effects, and, for binary outcomes, causal risk and odds ratios. Average treatment effects on treated and control units are additionally through the `att` and `atc` arguments.

# Examples of Missing Data

To illustrate the breadth of missingness patterns `drcmd` can accommodate, we outline a non-exhaustive set of example problems, assuming throughout that the causal assumptions of consistency, positivity, and unconfoundedness hold.

**Two-phase sampling.** The theory `drcmd` leverages is most often invoked in two-phase designs, where the researcher controls the missingness rate of difficult-to-measure variables. One observes $(Y_i, A_i, V_i, R_i U_i, R_i)$, where $R$ is a complete-case indicator and a subset $U$ of the confounders $X=(V,U)$ is measured only when $R_i=1$, with known selection probabilities $P(R=1\mid Y,A,V)$ that `drcmd` can incorporate directly.

**Missing outcome and treatment.** More broadly, one may observe $(R_i Y_i, R_i A_i, X_i, R_i)$ with $(Y,A)\perp R \mid X$; unlike two-phase designs, the mechanism $P(R=1\mid X)$ must be estimated.

**Surrogate outcomes.** One may observe $(R_i Y_i, A_i, X_i, R_i, S_i)$, where $S_i$ is an always-observed surrogate for the partially missing outcome $Y$; `drcmd` yields valid estimates when $Y \perp R \mid X, A, S$.

# Availability 
`drcmd` is publicly available for download on [GitHub](https://github.com/keithbarnatchez/drcmd). The package has been used to estimate causal effects under two-phase sampling with error-prone outcome and treatment measurements [@barnatchez2025twophase]. Further information on the use of the package can be found in the package vignette and user manual.

# AI usage disclosure

Large language models (Anthropic's Claude) were used to assist with debugging code and the construction of testing function.  All AI-assisted code was reviewed, verified, and edited by the authors, who take full responsibility for the the correctness of the software.

# References

