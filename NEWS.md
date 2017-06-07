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

# qte 1.2
  * Added ddid2 method to compute quantile treatment effects under a Difference in Differences assumption following the procedure of Callaway and Li, 2016.
  * Include covariates for all methods using formulas (xformla parameter) instead of passing in names
  * Improvements to bootstrapping standard errors