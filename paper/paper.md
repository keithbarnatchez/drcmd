---
title: "`drcmd`: Doubly-robust causal inference with missing data
tags:
    - causal inference
    - missing data
    - R
authors: 
    - name: Keith Barnatchez
affiliations:
    - name: Harvard T. H. Chan School of Public Health
    - index: 1
date: 16 May 2025
bibliography: refs.bib
---

# Summary

The `drcmd` `R` package performs semi-parametric efficient estimation of causal effects of point exposures in settings where the data available to the researcher is subject to missingness. By implementing methods from semi-parametric theory and missing data (see, e.g. @robins1994estimation and @kennedy2022semiparametric), `drcmd` accomodates general patterns of missing data, while enabling users to estimate nuisance functions with flexible machine leaening methods. `drcmd` automatically determines the missingness patterns present in user-supplied data, and provides information on assumptions that must hold regarding the missingness mechanisms in order for point estimates and inferences to be valid. By accommodating *general* patterns of missingness, `drcmd` serves as a centralized library for researchers aiming to perform causal inference with missing data.

# Statement of Need
The use of doubly-robust methods for performing causal inference of point exposures on outcomes of interest has surged over the past decade, and numerous software packages have been developed for implementing these estimators. While these packages are well-suited for use on complete, missingness-free data, leading statistical software packages provide little to no support for missing data. The lack of a centralized package for performing doubly-robust causal inference has functioned as a severe impediment for researchers, as missing data is ubiquiotous in real-world data, and the specific patterns of missingness can greatly vary across applications. `drcmd` addresses this shortcoming by providing a single R package for performing doubly-robust causal inference in the presence of missing data.

# Examples of Missing Data

### Example 1: Two phase sampling 

### Example 2: Surrogate outcomes

### Example 3: Missing treatments

<!--   -->
<!-- While a subset of existing -->
<!-- packages can address specific forms of missingness, there is currently no package that
accomodates estimation in the presence of missingness in (i) , (ii) and (iii). Given the
ubiquity of missing data in public health applications, there remains a critical need for
statistical software which can flexibly address missing data in a manner which lever-
ages modern developments at the intersection of causal inference and semiparametric
theory -->


# Availability 
`drcmd` is publically available for download on [GitHub](https://github.com/keithbarnatchez/drcmd). Further information on the use of the package can be found on --

# References

