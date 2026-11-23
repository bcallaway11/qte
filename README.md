

<!-- README.md is generated from README.qmd. Please edit that file. -->

# qte — Quantile Treatment Effects in R <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/qte.png)](https://CRAN.R-project.org/package=qte)
[![R-CMD-check](https://github.com/bcallaway11/qte/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/bcallaway11/qte/actions/workflows/R-CMD-check.yml)
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
- `panel_qtt()` — Panel QTT via copula stability (Callaway and Li 2019)
- `ddid()` — Distributional Difference-in-Differences (Callaway and Li
  2019)
- `mdid()` — Mean Difference-in-Differences (Thuysbaert 2007)
- `cdf_did()` — Difference-in-Differences in distribution functions
  (Roth and Sant’Anna 2023)
- `lou_qtt()` — Lagged-outcome unconfoundedness QTT

## Installation

``` r
# Install from CRAN:
install.packages("qte")

# Install the development version from GitHub:
# install.packages("remotes")
remotes::install_github("bcallaway11/qte")
```

## Quick start — unconfoundedness

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
  xformula   = xf,
  est_method = "aipw",
  target     = "qtt",
  probs      = seq(0.1, 0.9, 0.1),
  biters     = 100
)
summary(res_cs)
#> 
#> Overall ATT:  
#>        ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -4543.928      905.7607  -6319.186    -2768.67 *
#> 
#> 
#> QTT:
#>  Tau         QTT Std. Error [ 95% Simult.  Conf. Band]  
#>  0.1      0.0001   101.1285      -286.7827    286.7829  
#>  0.2  -1002.7420   835.9486     -3373.3467   1367.8627  
#>  0.3  -3400.6253  1684.5367     -8177.6779   1376.4273  
#>  0.4  -5009.2481  1254.5540     -8566.9448  -1451.5514 *
#>  0.5  -4602.5048  1025.3409     -7510.1929  -1694.8167 *
#>  0.6  -5229.1200  1295.7018     -8903.5046  -1554.7354 *
#>  0.7  -5505.9404  1261.1849     -9082.4409  -1929.4398 *
#>  0.8  -6885.7327  1383.7991    -10809.9459  -2961.5196 *
#>  0.9 -10517.0895  2558.9168    -17773.7318  -3260.4473 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Note: the outcome looks discrete (740 distinct values; 12% at 0); quantile effects may be unreliable.
```

Plot the QTT curve with a uniform confidence band:

``` r
autoplot(res_cs)
```

<img src="man/figures/README-unc-qte-plot-1.png" style="width:70.0%"
alt="QTT estimates across quantiles with uniform confidence band" />

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
#>  -0.0197        0.0151    -0.0495      0.0102 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0508     0.0213        0.0090      0.0927 *
#>          -2   0.0158     0.0152       -0.0139      0.0456  
#>          -1  -0.0128     0.0153       -0.0428      0.0172  
#>           0  -0.0081     0.0147       -0.0368      0.0207  
#>           1  -0.0364     0.0225       -0.0805      0.0077  
#>           2  -0.1226     0.0433       -0.2074     -0.0378 *
#>           3  -0.0930     0.0501       -0.1912      0.0053  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

Event-study plot showing pre-trends and post-treatment ATT by event
time:

``` r
autoplot(res_att, type = "dynamic")
```

<img src="man/figures/README-cic-att-plot-1.png" style="width:70.0%"
alt="Event-study plot of ATT estimates by event time" />

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
alt="Overall QTT curve from cic()" />

## Available estimators

| Function      | Method                             | Target     | Panel required |
|---------------|------------------------------------|------------|----------------|
| `unc_qte()`   | Unconfoundedness (IPW / OR / AIPW) | QTE or QTT | No             |
| `cic()`       | Change in Changes                  | ATT or QTT | Optional       |
| `qdid()`      | Quantile DiD                       | ATT or QTT | Optional       |
| `panel_qtt()` | Panel QTT (copula stability)       | ATT or QTT | Yes            |
| `ddid()`      | Distributional DiD                 | ATT or QTT | Yes            |
| `mdid()`      | Mean DiD                           | ATT or QTT | Optional       |
| `cdf_did()`   | DiD in distribution functions      | ATT or QTT | Optional       |
| `lou_qtt()`   | Lagged-outcome unconfoundedness    | ATT or QTT | Optional\*     |

All panel estimators support staggered treatment adoption and return
group-specific, event-study, and overall aggregations. “Optional” means
the estimator also works with repeated cross sections (`panel = FALSE`).
`ddid()` and `panel_qtt()` need the same units over time (a comparison
unit’s change, or the copula between periods). \*`lou_qtt()` needs panel
data for its default lagged-outcome covariate; with
`lagged_outcome_cov = FALSE` it works with repeated cross sections.

## Documentation and vignettes

Full documentation and vignettes are available at the [pkgdown
site](https://bcallaway11.github.io/qte/):

- **Quantile Treatment Effects in R** — `unc_qte()` under random
  assignment and selection on observables
- **Panel Data Estimators for Quantile Treatment Effects** —
  identification assumptions and usage for the panel estimators
- **Staggered Treatment Adoption** — applied workflow with `mpdta`: QTT
  curves, event-study plots, and cross-estimator comparison
