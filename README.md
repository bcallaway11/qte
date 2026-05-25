

<!-- README.md is generated from README.qmd. Please edit that file. -->

# qte — Quantile Treatment Effects in R <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/qte.png)](https://CRAN.R-project.org/package=qte)
[![R-CMD-check](https://github.com/bcallaway11/qte/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bcallaway11/qte/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The `qte` package provides methods for estimating Quantile Treatment
Effects (QTE) and Quantile Treatment Effects on the Treated (QTT) in R.
Where the average treatment effect summarizes the impact of a policy by
a single number, the QTE describes how treatment effects vary across the
outcome distribution — useful whenever the policy’s impact is
heterogeneous or when distributional consequences (e.g., for inequality)
are of interest.

**Cross-sectional estimators** (no panel data required):

- `unc_qte()` — QTE/QTT under unconfoundedness (IPW, outcome regression,
  or doubly robust); covers random assignment as a special case

**Panel and repeated cross-section estimators** (staggered treatment
adoption supported for all):

- `cic()` — Change in Changes (Athey and Imbens 2006)
- `qdid()` — Quantile Difference-in-Differences (Athey and Imbens 2006;
  Meyer, Viscusi, and Durbin 1995)
- `ddid()` — Distributional Difference-in-Differences (Callaway and Li
  2019)
- `mdid()` — Mean Difference-in-Differences (Thuysbaert 2007)
- `panel_qtt()` — Panel QTT via copula stability (Callaway and Li 2019)
- `lou_qtt()` — Lagged-outcome unconfoundedness QTT

## Installation

``` r
# Install from CRAN:
install.packages("qte")

# Install the development version from GitHub:
# install.packages("remotes")
remotes::install_github("bcallaway11/qte")
```

## Quick start — selection on observables

The `unc_qte()` function estimates the QTE or QTT under an
unconfoundedness assumption. Here we use the observational Lalonde
(1986) data to estimate the QTT of a job training program, controlling
for pre-treatment characteristics via doubly robust estimation.

``` r
data(lalonde)

xf <- ~ age + I(age^2) + education + black + hispanic + married + nodegree

res_cs <- unc_qte(
  yname      = "re78",
  dname      = "treat",
  data       = lalonde.psid,
  xformla    = xf,
  est_method = "aipw",
  target     = "qtt",
  probs      = seq(0.1, 0.9, 0.1),
  biters     = 100
)
summary(res_cs)
#> 
#> Overall ATT:  
#>        ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -4685.583       868.226  -6387.275   -2983.891 *
#> 
#> 
#> QTT:
#>  Tau         QTT Std. Error [ 95% Simult.  Conf. Band]  
#>  0.1      0.0001   152.1425      -298.1938    298.1940  
#>  0.2  -1002.7420   823.2062     -2616.1965    610.7125  
#>  0.3  -3400.5673  1666.5363     -6666.9184   -134.2163 *
#>  0.4  -5009.2491  1358.9281     -7672.6991  -2345.7990 *
#>  0.5  -4602.4652  1021.2978     -6604.1722  -2600.7583 *
#>  0.6  -5229.1454  1324.7191     -7825.5472  -2632.7437 *
#>  0.7  -5507.4720  1381.3583     -8214.8844  -2800.0595 *
#>  0.8  -6885.7529  1526.9716     -9878.5622  -3892.9436 *
#>  0.9 -10517.0625  2280.6228    -14987.0011  -6047.1240 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

Plot the QTT curve with a uniform confidence band:

``` r
autoplot(res_cs)
```

<img src="man/figures/README-unc-qte-plot-1.png" style="width:70.0%"
data-fig-alt="QTT estimates across quantiles with uniform confidence band" />

## Staggered treatment adoption

All panel estimators use a common `yname`/`gname`/`tname`/`idname`
interface and support staggered treatment adoption via
[ptetools](https://github.com/bcallaway11/ptetools). The example below
uses the `mpdta` dataset (county-level employment, from the `did`
package) with the Change in Changes estimator.

``` r
data(mpdta, package = "did")

res_att <- cic(
  yname   = "lemp",
  gname   = "first.treat",
  tname   = "year",
  idname  = "countyreal",
  data    = mpdta,
  gt_type = "att",
  biters  = 100
)
summary(res_att)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.] 
#>  -0.0197        0.0151     -0.057      0.0177 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0508     0.0213        0.0091      0.0926 *
#>          -2   0.0158     0.0156       -0.0148      0.0465  
#>          -1  -0.0128     0.0180       -0.0481      0.0224  
#>           0  -0.0081     0.0151       -0.0377      0.0216  
#>           1  -0.0364     0.0212       -0.0779      0.0051  
#>           2  -0.1226     0.0427       -0.2063     -0.0388 *
#>           3  -0.0930     0.0473       -0.1856     -0.0003 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

Event-study plot showing pre-trends and post-treatment ATT by event
time:

``` r
autoplot(res_att, type = "dynamic")
```

<img src="man/figures/README-cic-att-plot-1.png" style="width:70.0%"
data-fig-alt="Event-study plot of ATT estimates by event time" />

The same estimator returns a full QTT curve when `gt_type = "qtt"`:

``` r
res_qtt <- cic(
  yname   = "lemp",
  gname   = "first.treat",
  tname   = "year",
  idname  = "countyreal",
  data    = mpdta,
  gt_type = "qtt",
  probs   = seq(0.1, 0.9, 0.1),
  biters  = 100
)
autoplot(res_qtt)
```

<img src="man/figures/README-cic-qtt-1.png" style="width:70.0%"
data-fig-alt="Overall QTT curve from cic()" />

## Available estimators

| Function      | Method                             | Target     | Panel required |
|---------------|------------------------------------|------------|----------------|
| `unc_qte()`   | Unconfoundedness (IPW / OR / AIPW) | QTE or QTT | No             |
| `cic()`       | Change in Changes                  | ATT or QTT | Optional       |
| `qdid()`      | Quantile DiD                       | ATT or QTT | Optional       |
| `ddid()`      | Distributional DiD                 | ATT or QTT | Yes            |
| `mdid()`      | Mean DiD                           | ATT or QTT | Optional       |
| `panel_qtt()` | Panel QTT (copula stability)       | QTT        | Yes            |
| `lou_qtt()`   | Lagged-outcome unconfoundedness    | ATT or QTT | Yes            |

All panel estimators support staggered treatment adoption and return
group-specific, event-study, and overall aggregations.

## Documentation and vignettes

Full documentation and vignettes are available at the [pkgdown
site](https://bcallaway11.github.io/qte/):

- **Quantile Treatment Effects in R** — `unc_qte()` under random
  assignment and selection on observables
- **Panel Data Estimators for Quantile Treatment Effects** —
  identification assumptions and usage for all six panel estimators
- **Staggered Treatment Adoption** — applied workflow with `mpdta`: QTT
  curves, event-study plots, and cross-estimator comparison
