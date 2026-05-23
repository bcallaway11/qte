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
  weightsname = NULL,
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
  sections.

- xformula:

  One-sided formula for covariates used in the covariate adjustment.
  Default `~1` uses no covariates.

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
#>  -0.0197        0.0112    -0.0462      0.0069 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0508     0.0253        0.0012      0.1005 *
#>          -2   0.0158     0.0115       -0.0066      0.0383  
#>          -1  -0.0128     0.0172       -0.0466      0.0209  
#>           0  -0.0081     0.0104       -0.0285      0.0124  
#>           1  -0.0364     0.0214       -0.0784      0.0056  
#>           2  -0.1226     0.0436       -0.2080     -0.0371 *
#>           3  -0.0930     0.0412       -0.1736     -0.0123 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## Full QTT curve at selected quantiles
res_qtt <- cic(yname = "lemp", gname = "first.treat", tname = "year",
               idname = "countyreal", data = mpdta,
               gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res_qtt)
#> 
#> Overall QTT Curve:
#>  Quantile     QTT Std. Error 95% CB Lower 95% CB Upper
#>       0.1  0.0286     0.1817      -0.5410       0.5982
#>       0.2 -0.0469     0.0660      -0.2536       0.1599
#>       0.3 -0.0518     0.0641      -0.2527       0.1490
#>       0.4  0.0127     0.0498      -0.1434       0.1688
#>       0.5 -0.0455     0.0541      -0.2149       0.1240
#>       0.6 -0.0445     0.0398      -0.1692       0.0802
#>       0.7  0.0046     0.0497      -0.1512       0.1604
#>       0.8 -0.0187     0.0537      -0.1871       0.1497
#>       0.9  0.0166     0.0681      -0.1968       0.2300
#> 
# }
```
