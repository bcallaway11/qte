---
output: github_document
title: "NEWS"
author: "Brantly Callaway"
date: "`r Sys.Date()`"
vignette: >
  %\VignetteIndexEntry{NEWS}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# qte 1.5.0

  * Renamed `ci.qte()` to `unc_qte()` with a modernised argument API:
    `yname`/`dname` replace the formula interface, `weightsname` replaces
    the `w` vector, `biters`/`cl`/`alp`/`cband` match the `did`/`ptetools`
    conventions. `ci.qte()` is retained as a deprecated wrapper.
  * New `cic()` function: unified interface for the Change in Changes
    estimator covering two-period and staggered treatment adoption, for
    both panel and repeated cross sections data.  Replaces the older
    `CiC()` (now deprecated) and `cic2()` (removed).
  * New `cic_gt()`: the 2x2 CiC workhorse used as a `ptetools` attgt_fun.
    Automatically detects panel vs. repeated cross sections from the data
    and returns counterfactual distributions for downstream aggregation.
  * `cic()` now supports `gt_type = "qtt"` to return a full Quantile
    Treatment Effect on the Treated curve at each aggregation level
    (overall, group-specific, event-study).  Accepts a `probs` argument
    to control the quantile grid.  CDFs are mixed across group-time cells
    before inversion, avoiding bias from averaging scalar quantiles.
  * Removed `qtes2mat()` and `diffQ()`: formatting/utility helpers tied to
    the legacy `QTE` class with no internal callers. No replacement planned.
  * Removed `DR` (doubly robust ATT estimator, unused).
  * Removed `spatt` (Abadie 2005 semiparametric DiD ATT estimator).
    Users needing a semiparametric ATT with staggered adoption support
    should use `did::att_gt()` instead.
  * Removed `bounds` (Fan-Yu 2012 partial identification bounds for the
    QTT). The function only supported two periods, did not correctly handle
    covariates, and had no sampling weights support. No replacement is
    planned.
  * Added smoke tests for `cic()` covering panel ATT, repeated cross
    sections ATT, and panel QTT.
  * Added vignette demonstrating `cic()` with staggered treatment adoption
    using the `mpdta` dataset.

# qte 1.4.0
  * Added new function `cic2` for change-in-changes with multiple
    periods and variation in treatment timing.  This uses the `pte`
	package as a backend, and is not fully tested yet.

# qte 1.3.3
  * bug fixes to allow critical values for uniform confidence 
    band to be affected by `alp` argument rather than stuck 
	on `alp=0.05`

# qte 1.3.2
  * now pass the significance level through in the results
  * added functionality to change significance level in `ggqte` in 
    for pointwise confidence intervals w/o needing to re-estimate

# qte 1.3.1
  * Updated maintainer's contact information
  * Improvements/updates to vignettes
  
# qte 1.3.0
  * Added first-step quantile regression support for panel.qtet, QDiD, CiC
  * Added some checks that inputted data is of correct format, particularly in panel data case
  
# qte 1.2.3 
  * Now using empirical quantiles in CiC method
  * Added code for computing uniform confidence bands
  * General improvements to codebase, e.g., combining same code used across methods

# qte 1.2.2 
  * Some bug fixes to CiC method

# qte 1.2
  * Added ddid2 method to compute quantile treatment effects under a Difference in Differences assumption following the procedure of Callaway and Li, 2016.
  * Include covariates for all methods using formulas (xformla parameter) instead of passing in names
  * Improvements to bootstrapping standard errors
  * Some extensions to computing the ATT by itself in \code{spatt}
  * Preliminary version of multi period ATT moved to \code{did} package
