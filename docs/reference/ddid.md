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
  weightsname = NULL,
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

- weightsname:

  Name of the column in `data` containing sampling weights. Default
  `NULL` uses equal weights.

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
#>  -0.0453        0.0137     -0.072     -0.0186 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0134     0.0145       -0.0151      0.0419  
#>          -2  -0.0146     0.0122       -0.0385      0.0094  
#>          -1  -0.0403     0.0147       -0.0692     -0.0114 *
#>           0  -0.0324     0.0141       -0.0600     -0.0048 *
#>           1  -0.0637     0.0168       -0.0967     -0.0308 *
#>           2  -0.1377     0.0418       -0.2197     -0.0557 *
#>           3  -0.1087     0.0345       -0.1763     -0.0410 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## Full QTT curve at selected quantiles
res_qtt <- ddid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res_qtt)
#> 
#> Overall QTT Curve:
#>  Quantile     QTT Std. Error 95% CB Lower 95% CB Upper
#>       0.1  0.0286     0.1354      -0.7916       0.8487
#>       0.2 -0.0620     0.0536      -0.3867       0.2627
#>       0.3 -0.0624     0.0426      -0.3206       0.1957
#>       0.4 -0.0601     0.0364      -0.2806       0.1604
#>       0.5 -0.0614     0.0383      -0.2931       0.1704
#>       0.6 -0.0190     0.0384      -0.2518       0.2137
#>       0.7 -0.0064     0.0473      -0.2929       0.2802
#>       0.8 -0.0540     0.0449      -0.3257       0.2177
#>       0.9 -0.0185     0.0763      -0.4803       0.4433
#> 
# }
```
