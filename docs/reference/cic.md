# Change in Changes

Computes Quantile Treatment effects on the Treated (QTT) and the Average
Treatment Effect on the Treated (ATT) using the Change in Changes
identification strategy of Athey and Imbens (2006). Handles two-period
data and staggered treatment adoption uniformly: a two-period, two-group
dataset is the degenerate single-(g,t) case. Supports both panel and
repeated cross sections data.

## Usage

``` r
cic(
  yname,
  gname,
  tname,
  idname = NULL,
  data,
  panel = TRUE,
  xformula = ~1,
  cond_dist = c("qr", "dr"),
  weightsname = NULL,
  clustervars = NULL,
  control_group = "notyettreated",
  anticipation = 0,
  alp = 0.05,
  cband = TRUE,
  biters = 100,
  cl = 1,
  ret_quantile = NULL,
  gt_type = "att",
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

- panel:

  Logical; `TRUE` (default) for panel data, `FALSE` for repeated cross
  sections (see the section on repeated cross sections).

- xformula:

  One-sided formula for covariates used in the covariate adjustment.
  Default `~1` uses no covariates. With panel data, each unit's
  covariates are taken from the base period of each group-time cell
  (never from a post-treatment period); with repeated cross sections,
  each observation's own values are used.

- cond_dist:

  How the conditional distributions are estimated when `xformula` has
  covariates: `"qr"` (default) for quantile regression at the
  percentiles 0.01, ..., 0.99, or `"dr"` for distribution regression (a
  logit of \\1\\Y \le y\\\\ on the covariates at a grid of thresholds;
  Chernozhukov, Fernandez-Val and Melly 2013). Distribution regression
  handles outcomes with mass points (e.g. many zeros) naturally. Ignored
  without covariates.

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

  Significance level for confidence bands. Default `0.05`.

- cband:

  Logical; if `TRUE` (default) compute a uniform confidence band rather
  than pointwise intervals.

- biters:

  Number of bootstrap iterations. Default `100`.

- cl:

  Number of clusters for parallel computation. Default `1`.

- ret_quantile:

  Passed through to `ptetools` for the `"qott"` case. Ignored when
  `gt_type = "qtt"` (use `probs` instead).

- gt_type:

  Type of group-time effect to compute. `"att"` (default) returns
  ATT(g,t). `"qtt"` returns the full QTT curve over `probs` using
  mixture-CDF aggregation. `"qott"` returns the quantile of the
  individual treatment effect distribution under rank invariance (panel
  only).

- probs:

  For `gt_type = "qtt"`, the quantile grid at which to evaluate the QTT
  curve. Default is `seq(0.05, 0.95, 0.05)`.

## Value

For `gt_type = "att"`, a `pte_results` object from `ptetools`. For
`gt_type = "qtt"`, a `pte_qtt` object with overall, group-specific, and
dynamic QTT curves and bootstrap SEs.

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

`cic()`,
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
(the copula between periods), the lagged outcome in
[`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md),
and `gt_type = "qott"` (the distribution of individual treatment
effects).

## References

Athey, Susan and Guido Imbens. “Identification and Inference in
Nonlinear Difference-in-Differences Models.” Econometrica 74(2), pp.
431-497, 2006.

## Examples

``` r
# \donttest{
data(mpdta, package = "did")

## ATT aggregated across all groups and periods
res_att <- cic(yname = "lemp", gname = "first.treat", tname = "year",
               idname = "countyreal", data = mpdta,
               gt_type = "att", biters = 20)
summary(res_att)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.] 
#>  -0.0197        0.0204    -0.0565      0.0172 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0508     0.0189        0.0137      0.0880 *
#>          -2   0.0158     0.0153       -0.0142      0.0459  
#>          -1  -0.0128     0.0188       -0.0497      0.0240  
#>           0  -0.0081     0.0197       -0.0467      0.0306  
#>           1  -0.0364     0.0234       -0.0823      0.0095  
#>           2  -0.1226     0.0419       -0.2047     -0.0404 *
#>           3  -0.0930     0.0377       -0.1669     -0.0190 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## Full QTT curve at selected quantiles
res_qtt <- cic(yname = "lemp", gname = "first.treat", tname = "year",
               idname = "countyreal", data = mpdta,
               gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res_qtt)
#> 
#> Overall QTT curve:
#>  Quantile     QTT Std. Error [95% Simult. Conf. Band]        
#>       0.1  0.0278     0.1125                  -0.2834 0.3389 
#>       0.2 -0.0230     0.0871                  -0.2638 0.2178 
#>       0.3 -0.0103     0.0630                  -0.1845 0.1639 
#>       0.4 -0.0690     0.0553                  -0.2218 0.0838 
#>       0.5 -0.0185     0.0653                  -0.1990 0.1620 
#>       0.6 -0.0149     0.0373                  -0.1180 0.0882 
#>       0.7  0.0389     0.0486                  -0.0955 0.1732 
#>       0.8 -0.0387     0.0476                  -0.1702 0.0929 
#>       0.9  0.0135     0.0736                  -0.1901 0.2170 
#> 
#> Dynamic effects (quantile 0.5):
#>  Event Time     QTT Std. Error [95% Simult. Conf. Band]        
#>          -3  0.0159     0.0528                  -0.1506 0.1825 
#>          -2  0.0868     0.0774                  -0.1847 0.3582 
#>          -1 -0.0683     0.0762                  -0.2474 0.1109 
#>           0 -0.0237     0.0577                  -0.2029 0.1555 
#>           1 -0.0429     0.0512                  -0.1774 0.0915 
#>           2 -0.0960     0.1165                  -0.4172 0.2251 
#>           3 -0.1972     0.1554                  -0.6255 0.2312 
#> 
#> Group effects (quantile 0.5):
#>  Group     QTT Std. Error [95% Simult. Conf. Band]        
#>   2004 -0.1391     0.1101                  -0.4999 0.2218 
#>   2006  0.0149     0.0673                  -0.1872 0.2170 
#>   2007 -0.1181     0.0557                  -0.2993 0.0630 
#> 
#> ATT, and test of QTT(tau) = ATT for all quantiles:
#>      Parameter     ATT Std. Error [95% Pointwise Conf. Int.] p (QTT = ATT)
#>        Overall -0.0197     0.0164        -0.0518      0.0124          0.95
#>  Event time -3  0.0508     0.0229         0.0059      0.0957          0.95
#>  Event time -2  0.0158     0.0138        -0.0112      0.0429          0.90
#>  Event time -1 -0.0128     0.0166        -0.0454      0.0198          0.95
#>   Event time 0 -0.0081     0.0175        -0.0423      0.0262          0.85
#>   Event time 1 -0.0364     0.0226        -0.0807      0.0079          0.90
#>   Event time 2 -0.1226     0.0482        -0.2169     -0.0282          0.85
#>   Event time 3 -0.0930     0.0469        -0.1849     -0.0010          0.95
#> 
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Simultaneous bands are uniform over quantiles, for each QTT curve.
#> p (QTT = ATT): bootstrap sup-t test, from the same draws, that QTT(tau) - ATT = 0
#> at every estimated quantile.
# }
```
