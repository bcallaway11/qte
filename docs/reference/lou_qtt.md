# Lagged Outcome Unconfoundedness QTT

Estimates the Quantile Treatment Effect on the Treated (QTT) and Average
Treatment Effect on the Treated (ATT) under a lagged-outcome
unconfoundedness assumption with staggered treatment adoption. The key
identifying assumption is \\Y\_{g,t}(0) \perp D \mid X,
Y\_{\text{pre}}\\, i.e., conditional on observed covariates and the
pre-treatment outcome, treatment is as good as randomly assigned within
each cohort-period cell.

Estimation operates at the (g,t) level: for each cohort \\g\\ and
post-treatment period \\t\\, a cross-sectional comparison is made
between the treated group and a not-yet-treated (or never-treated)
comparison group, adjusting for covariates and (optionally) the
pre-treatment outcome. Group-time estimates are then aggregated to
overall, dynamic, and group-specific summaries.

Three estimation methods are available:

- `"ipw"`:

  Propensity-score reweighting. Control units are reweighted by \\\hat
  p(X,Y\_{\text{pre}})/(1-\hat p(X,Y\_{\text{pre}}))\\ to approximate
  the covariate distribution of the treated group.

- `"reg"`:

  Outcome regression. A quantile regression model is fit on control
  units' post-period outcomes as a function of \\(X, Y\_{\text{pre}})\\,
  then predicted at treated units to construct the counterfactual
  distribution (Melly 2006; Chernozhukov, Fernandez-Val, and Melly
  2013).

- `"aipw"`:

  Doubly-robust augmented IPW. Combines the propensity score and outcome
  models. Consistent if either model is correctly specified.

When `lagged_outcome_cov = FALSE` and `xformula = ~1`, all three methods
reduce to a simple distribution comparison within each (g,t) cell, which
is consistent under unconditional unconfoundedness.

The ATT analogue (without the QTT) for staggered adoption under
lagged-outcome unconfoundedness is developed in Callaway (2023).

## Usage

``` r
lou_qtt(
  yname,
  gname,
  tname,
  idname = NULL,
  data,
  xformula = ~1,
  lagged_outcome_cov = TRUE,
  est_method = c("ipw", "reg", "aipw"),
  cond_dist = c("qr", "dr"),
  panel = TRUE,
  weightsname = NULL,
  clustervars = NULL,
  control_group = "notyettreated",
  anticipation = 0,
  alp = 0.05,
  cband = TRUE,
  biters = 100,
  cl = 1,
  gt_type = c("att", "qtt"),
  probs = NULL
)
```

## Arguments

- yname:

  Name of the outcome variable in `data`.

- gname:

  Name of the treatment group variable (first treatment period; 0 for
  never-treated units).

- tname:

  Name of the time period variable.

- idname:

  Name of the unit id variable. Required when `panel = TRUE`.

- data:

  A data frame.

- xformula:

  One-sided formula for additional covariates used in the propensity
  score and/or outcome model. Default `~1` uses no additional covariates
  beyond the lagged outcome (when `lagged_outcome_cov = TRUE`). With
  panel data, each unit's covariates are taken from the base period of
  each group-time cell (never from a post-treatment period); with
  repeated cross sections, each observation's own values are used.

- lagged_outcome_cov:

  Logical; if `TRUE` (default), the pre-treatment outcome is added as a
  covariate in both the propensity score and outcome models. Requires
  `panel = TRUE`.

- est_method:

  Estimation method: `"ipw"` (default), `"reg"` (`"or"` is a deprecated
  alias), or `"aipw"`. See Details.

- cond_dist:

  For `est_method = "reg"` or `"aipw"`, how the conditional distribution
  of the untreated outcome is estimated: `"qr"` (default) for quantile
  regression or `"dr"` for distribution regression. See
  [`cic_gt`](https://bcallaway11.github.io/qte/reference/cic_gt.md).
  Ignored for `"ipw"`.

- panel:

  Logical; `TRUE` (default) for panel data, `FALSE` for repeated cross
  sections. `lagged_outcome_cov = TRUE` requires `panel = TRUE`.

- weightsname:

  Name of the column in `data` containing sampling weights. Default
  `NULL` uses equal weights.

- clustervars:

  Name of a column in `data` to cluster standard errors on, in addition
  to `idname` (as in
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)):
  one variable, constant over time within each unit. Default `NULL`
  clusters at the unit level. With the empirical bootstrap, whole
  clusters are resampled.

- control_group:

  Which units to use as the comparison group: `"notyettreated"`
  (default) or `"nevertreated"`.

- anticipation:

  Number of periods of anticipation. Default `0`.

- alp:

  Significance level for confidence intervals. Default `0.05`.

- cband:

  Logical; if `TRUE` (default), report a simultaneous confidence band
  (uniform over all quantiles in `probs`) in addition to pointwise
  intervals.

- biters:

  Number of bootstrap iterations. Default `100`.

- cl:

  Number of clusters for parallel bootstrap. Default `1`.

- gt_type:

  Type of group-time effect to compute. `"att"` (default) returns the
  ATT curve aggregated over (g,t) cells. `"qtt"` returns the full QTT
  curve over `probs` using mixture-CDF aggregation.

- probs:

  For `gt_type = "qtt"`, the quantile grid at which to evaluate the QTT.
  Default `seq(0.05, 0.95, 0.05)`.

## Value

For `gt_type = "att"`, a `pte_emp_boot` object. For `gt_type = "qtt"`, a
`pte_qtt` object with overall, group-specific, and dynamic QTT curves,
bootstrap standard errors, and pointwise and uniform confidence bands.

## Discrete outcomes

Quantile treatment effects usually require a continuous outcome. With
mass points (e.g. many zeros) or few distinct values, quantiles jump and
effects at quantiles inside a mass point are 0 by construction;
identification based on ranks or quantile functions (change-in-changes,
quantile DiD, the panel and distributional DiD copula arguments) only
bounds the effects (Athey and Imbens 2006, Section 4); and bootstrap
inference for quantiles can be unreliable. The estimators warn (once per
call, class `qte_discrete_outcome`) when more than 10\\ outcome's
observations share one value or it takes at most 20 distinct values,
across all groups and periods; the thresholds are the options
`qte.mass_point_tol` and `qte.discrete_n_unique`. The diagnostics are
returned as `outcome_discreteness`.
[`cdf_did`](https://bcallaway11.github.io/qte/reference/cdf_did.md) is
identified for discrete outcomes.

## Repeated cross sections

[`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
[`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
[`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md),
[`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
and `lou_qtt(lagged_outcome_cov = FALSE)` also work with repeated cross
sections (`panel = FALSE`; `idname` not needed). Without covariates they
use only each period's outcome distributions, so on panel data the
estimates are the same either way; the standard errors differ, because
the bootstrap then resamples observations rather than units. Some things
need the same units in several periods and so require panel data:
[`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) (each
comparison unit's change over time),
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
(the copula between periods), the lagged outcome in `lou_qtt()`, and
`gt_type = "qott"` (the distribution of individual treatment effects).

## References

Callaway, Brantly. “Policy Evaluation during a Pandemic.” *Journal of
Econometrics* 236(2), 2023.

Melly, Blaise. “Estimation of Counterfactual Distributions Using
Quantile Regression.” Working paper, University of St. Gallen, 2006.

Chernozhukov, Victor, Ivan Fernandez-Val, and Blaise Melly. “Inference
on Counterfactual Distributions.” *Econometrica* 81(6), pp. 2205–2268,
2013.

## See also

[`unc_qte`](https://bcallaway11.github.io/qte/reference/unc_qte.md) for
the cross-sectional (non-staggered) version.
[`cic`](https://bcallaway11.github.io/qte/reference/cic.md),
[`qdid`](https://bcallaway11.github.io/qte/reference/qdid.md),
[`mdid`](https://bcallaway11.github.io/qte/reference/mdid.md) for
alternative identification strategies with staggered adoption.

## Examples

``` r
# \donttest{
data(mpdta, package = "did")

## ATT under lagged-outcome unconfoundedness (IPW with pre-period outcome)
res_att <- lou_qtt(
  yname = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  lagged_outcome_cov = TRUE, est_method = "ipw",
  gt_type = "att", biters = 20
)
summary(res_att)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -0.0384        0.0177    -0.0727      -0.004 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0294     0.0209       -0.0115      0.0704  
#>          -2  -0.0124     0.0192       -0.0501      0.0252  
#>          -1  -0.0310     0.0198       -0.0699      0.0079  
#>           0  -0.0243     0.0168       -0.0573      0.0086  
#>           1  -0.0869     0.0337       -0.1530     -0.0208 *
#>           2  -0.1394     0.0344       -0.2069     -0.0719 *
#>           3  -0.1093     0.0412       -0.1901     -0.0284 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## QTT with doubly-robust estimation
res_qtt <- lou_qtt(
  yname = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  lagged_outcome_cov = TRUE, est_method = "aipw",
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20
)
summary(res_qtt)
#> 
#> Overall QTT curve:
#>  Quantile     QTT Std. Error [95% Simult. Conf. Band]        
#>       0.1  0.0564     0.1809                  -0.5888 0.7015 
#>       0.2 -0.0230     0.0530                  -0.2120 0.1660 
#>       0.3 -0.0103     0.0410                  -0.1563 0.1357 
#>       0.4 -0.0653     0.0410                  -0.2114 0.0809 
#>       0.5 -0.0237     0.0509                  -0.2052 0.1578 
#>       0.6 -0.0278     0.0419                  -0.1773 0.1218 
#>       0.7 -0.0025     0.0460                  -0.1665 0.1614 
#>       0.8 -0.0227     0.0556                  -0.2208 0.1754 
#>       0.9 -0.0613     0.0769                  -0.3356 0.2130 
#> 
#> Dynamic effects (quantile 0.5):
#>  Event Time     QTT Std. Error [95% Simult. Conf. Band]        
#>          -3  0.0192     0.0584                  -0.1552 0.1936 
#>          -2  0.0561     0.0619                  -0.1439 0.2561 
#>          -1 -0.0385     0.0398                  -0.1676 0.0906 
#>           0 -0.0237     0.0533                  -0.2201 0.1727 
#>           1 -0.0105     0.0815                  -0.2579 0.2370 
#>           2 -0.2735     0.1491                  -0.8034 0.2564 
#>           3 -0.3588     0.2177                  -1.1415 0.4239 
#> 
#> Group effects (quantile 0.5):
#>  Group     QTT Std. Error [95% Simult. Conf. Band]        
#>   2004 -0.2735     0.1094                  -0.6855 0.1386 
#>   2006 -0.0082     0.0747                  -0.2862 0.2698 
#>   2007 -0.0944     0.0609                  -0.2902 0.1013 
#> 
#> ATT, and test of QTT(tau) = ATT for all quantiles:
#>      Parameter     ATT Std. Error [95% Pointwise Conf. Int.] p (QTT = ATT)
#>        Overall -0.0297     0.0127        -0.0545     -0.0049          0.90
#>  Event time -3  0.0287     0.0137         0.0019      0.0556          1.00
#>  Event time -2 -0.0044     0.0138        -0.0313      0.0226          0.45
#>  Event time -1 -0.0222     0.0172        -0.0559      0.0114          0.45
#>   Event time 0 -0.0180     0.0120        -0.0414      0.0055          1.00
#>   Event time 1 -0.0523     0.0161        -0.0839     -0.0207          0.85
#>   Event time 2 -0.1380     0.0317        -0.2001     -0.0760          0.40
#>   Event time 3 -0.1066     0.0306        -0.1666     -0.0466          0.55
#> 
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Simultaneous bands are uniform over quantiles, for each QTT curve.
#> p (QTT = ATT): bootstrap sup-t test, from the same draws, that QTT(tau) - ATT = 0
#> at every estimated quantile.
# }
```
