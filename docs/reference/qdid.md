# Quantile Difference-in-Differences

Computes Quantile Treatment effects on the Treated (QTT) and the Average
Treatment Effect on the Treated (ATT) using the Quantile
Difference-in-Differences identification strategy of Athey and Imbens
(2006). Handles two-period data and staggered treatment adoption
uniformly. Supports both panel and repeated cross sections data.

## Usage

``` r
qdid(
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

  How conditional distributions are estimated when `xformula` has
  covariates: `"qr"` (default) for quantile regression or `"dr"` for
  distribution regression. See
  [`cic_gt`](https://bcallaway11.github.io/qte/reference/cic_gt.md).
  Ignored without covariates.

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

  Passed through to `ptetools` for the `"qott"` case.

- gt_type:

  Type of group-time effect to compute. `"att"` (default) returns
  ATT(g,t). `"qtt"` returns the full QTT curve over `probs` using
  mixture-CDF aggregation. `"qott"` returns the quantile of the
  individual treatment effect distribution (panel only).

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

[`cic()`](https://bcallaway11.github.io/qte/reference/cic.md), `qdid()`,
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
res_att <- qdid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "att", biters = 20)
summary(res_att)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.] 
#>  -0.0282        0.0137    -0.0567       3e-04 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0266     0.0187       -0.0100      0.0631  
#>          -2  -0.0040     0.0126       -0.0286      0.0206  
#>          -1  -0.0247     0.0138       -0.0516      0.0023  
#>           0  -0.0164     0.0137       -0.0434      0.0105  
#>           1  -0.0584     0.0199       -0.0974     -0.0194 *
#>           2  -0.1470     0.0328       -0.2113     -0.0827 *
#>           3  -0.1111     0.0317       -0.1734     -0.0489 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## Full QTT curve at selected quantiles
res_qtt <- qdid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res_qtt)
#> 
#> Overall QTT curve:
#>  Quantile     QTT Std. Error [95% Simult. Conf. Band]        
#>       0.1  0.0227     0.0739                  -0.1957 0.2411 
#>       0.2 -0.0348     0.0555                  -0.1989 0.1293 
#>       0.3 -0.0150     0.0484                  -0.1583 0.1282 
#>       0.4 -0.0807     0.0478                  -0.2221 0.0606 
#>       0.5 -0.0406     0.0356                  -0.1459 0.0646 
#>       0.6 -0.0080     0.0404                  -0.1273 0.1113 
#>       0.7  0.0325     0.0439                  -0.0973 0.1622 
#>       0.8  0.0010     0.0495                  -0.1453 0.1473 
#>       0.9 -0.0368     0.0394                  -0.1534 0.0797 
#> 
#> Dynamic effects (quantile 0.5):
#>  Event Time     QTT Std. Error [95% Simult. Conf. Band]        
#>          -3  0.0347     0.0675                  -0.1774 0.2468 
#>          -2  0.0523     0.0757                  -0.1962 0.3008 
#>          -1 -0.0023     0.0366                  -0.1166 0.1120 
#>           0 -0.0406     0.0581                  -0.2150 0.1337 
#>           1 -0.0679     0.0505                  -0.2258 0.0900 
#>           2 -0.0997     0.0918                  -0.4158 0.2164 
#>           3 -0.1691     0.0874                  -0.4643 0.1262 
#> 
#> Group effects (quantile 0.5):
#>  Group     QTT Std. Error [95% Simult. Conf. Band]        
#>   2004 -0.1110     0.0622                  -0.3404 0.1185 
#>   2006 -0.0324     0.0497                  -0.1988 0.1340 
#>   2007 -0.1319     0.0766                  -0.3650 0.1012 
#> 
#> ATT, and test of QTT(tau) = ATT for all quantiles:
#>      Parameter     ATT Std. Error [95% Pointwise Conf. Int.] p (QTT = ATT)
#>        Overall -0.0282     0.0097        -0.0472     -0.0093          0.95
#>  Event time -3  0.0266     0.0117         0.0037      0.0494          0.85
#>  Event time -2 -0.0040     0.0108        -0.0251      0.0171          0.80
#>  Event time -1 -0.0247     0.0141        -0.0522      0.0029          0.95
#>   Event time 0 -0.0164     0.0093        -0.0347      0.0019          0.70
#>   Event time 1 -0.0584     0.0241        -0.1057     -0.0111          0.95
#>   Event time 2 -0.1470     0.0440        -0.2333     -0.0607          0.20
#>   Event time 3 -0.1111     0.0478        -0.2048     -0.0175          0.05
#> 
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Simultaneous bands are uniform over quantiles, for each QTT curve.
#> p (QTT = ATT): bootstrap sup-t test, from the same draws, that QTT(tau) - ATT = 0
#> at every estimated quantile.
# }
```
