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
