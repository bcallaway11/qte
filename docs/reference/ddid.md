# Distributional Difference-in-Differences

Computes Quantile Treatment effects on the Treated (QTT) and the Average
Treatment Effect on the Treated (ATT) using the distributional DiD
identification strategy of Callaway, Li, and Oka (2018). Handles
two-period data and staggered treatment adoption uniformly via
`ptetools`. Requires panel data.

**Identification.** Under distributional parallel trends and a copula
restriction (the rank correlation of untreated potential outcomes
between the treated and control groups in the pre-period is preserved),
the counterfactual distribution \\F\_{Y(0),\text{post}\|D=1}\\ is
recovered by adding each control unit's actual change \\\Delta
Y\_{\text{ctrl},j}\\ to the treated pre-period quantile at that control
unit's rank. This point-identifies the QTT with only two time periods,
in contrast to methods (e.g. `panel_qtt`) that require three periods.

**Panel data required.** Unlike `cic`, `qdid`, and `mdid`, this
estimator cannot be applied to repeated cross sections because it
requires observing the within-unit outcome change for each control unit.

## Usage

``` r
ddid(
  yname,
  gname,
  tname,
  idname,
  data,
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

  Name of the unit id variable (required).

- data:

  A data frame.

- xformula:

  One-sided formula for covariates. Default `~1` uses no covariates.
  Each unit's covariates are taken from the base period of each
  group-time cell.

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

- gt_type:

  Type of group-time effect to compute. `"att"` (default) returns
  ATT(g,t). `"qtt"` returns the full QTT curve over `probs` using
  mixture-CDF aggregation.

- probs:

  For `gt_type = "qtt"`, the quantile grid. Default is
  `seq(0.05, 0.95, 0.05)`.

## Value

For `gt_type = "att"`, a `pte_results` object from `ptetools`. For
`gt_type = "qtt"`, a `pte_qtt` object.

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
`ddid()` (each comparison unit's change over time),
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
(the copula between periods), the lagged outcome in
[`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md),
and `gt_type = "qott"` (the distribution of individual treatment
effects).

## References

Callaway, Brantly, Tong Li, and Tatsushi Oka. “Quantile Treatment
Effects in Difference in Differences Models under Dependence
Restrictions and with Only Two Time Periods.” Journal of Econometrics
206(2), pp. 395-413, 2018.

## Examples

``` r
# \donttest{
data(mpdta, package = "did")

## ATT aggregated across all groups and periods
res_att <- ddid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "att", biters = 20)
summary(res_att)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -0.0453        0.0129    -0.0712     -0.0194 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0134     0.0175       -0.0208      0.0476  
#>          -2  -0.0146     0.0106       -0.0353      0.0061  
#>          -1  -0.0403     0.0124       -0.0646     -0.0160 *
#>           0  -0.0324     0.0118       -0.0556     -0.0092 *
#>           1  -0.0637     0.0174       -0.0978     -0.0296 *
#>           2  -0.1377     0.0483       -0.2324     -0.0431 *
#>           3  -0.1087     0.0369       -0.1810     -0.0363 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## Full QTT curve at selected quantiles
res_qtt <- ddid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res_qtt)
#> 
#> Overall QTT curve:
#>  Quantile     QTT Std. Error [95% Simult. Conf. Band]        
#>       0.1 -0.0062     0.1101                  -0.3219 0.3094 
#>       0.2 -0.0931     0.0623                  -0.2717 0.0854 
#>       0.3 -0.0156     0.0431                  -0.1392 0.1080 
#>       0.4 -0.0707     0.0423                  -0.1920 0.0506 
#>       0.5 -0.0884     0.0430                  -0.2117 0.0350 
#>       0.6 -0.0188     0.0383                  -0.1286 0.0910 
#>       0.7 -0.0033     0.0406                  -0.1196 0.1130 
#>       0.8 -0.0070     0.0605                  -0.1802 0.1663 
#>       0.9 -0.0661     0.0578                  -0.2318 0.0997 
#> 
#> Dynamic effects (quantile 0.5):
#>  Event Time     QTT Std. Error [95% Simult. Conf. Band]        
#>          -3  0.0201     0.0543                  -0.1842 0.2244 
#>          -2  0.0588     0.0457                  -0.0860 0.2037 
#>          -1 -0.0637     0.0680                  -0.2426 0.1152 
#>           0 -0.0880     0.0454                  -0.2443 0.0682 
#>           1 -0.0303     0.0503                  -0.2174 0.1568 
#>           2 -0.2837     0.2217                  -1.0261 0.4587 
#>           3 -0.3758     0.2062                  -1.1838 0.4323 
#> 
#> Group effects (quantile 0.5):
#>  Group     QTT Std. Error [95% Simult. Conf. Band]        
#>   2004 -0.2913     0.1805                  -0.7867 0.2040 
#>   2006 -0.0266     0.0517                  -0.2198 0.1667 
#>   2007 -0.1231     0.0537                  -0.2584 0.0122 
#> 
#> ATT, and test of QTT(tau) = ATT for all quantiles:
#>      Parameter     ATT Std. Error [95% Pointwise Conf. Int.] p (QTT = ATT)
#>        Overall -0.0453     0.0159        -0.0765     -0.0140          1.00
#>  Event time -3  0.0134     0.0137        -0.0134      0.0401          1.00
#>  Event time -2 -0.0146     0.0104        -0.0349      0.0058          0.80
#>  Event time -1 -0.0403     0.0158        -0.0713     -0.0094          1.00
#>   Event time 0 -0.0324     0.0157        -0.0631     -0.0016          0.95
#>   Event time 1 -0.0637     0.0182        -0.0993     -0.0281          0.90
#>   Event time 2 -0.1377     0.0365        -0.2093     -0.0662          0.85
#>   Event time 3 -0.1087     0.0444        -0.1958     -0.0215          0.60
#> 
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Simultaneous bands are uniform over quantiles, for each QTT curve.
#> p (QTT = ATT): bootstrap sup-t test, from the same draws, that QTT(tau) - ATT = 0
#> at every estimated quantile.
# }
```
