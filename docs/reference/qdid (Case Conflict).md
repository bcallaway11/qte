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
#>  -0.0271         0.016    -0.0581      0.0039 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0266     0.0195       -0.0117      0.0649  
#>          -2  -0.0027     0.0138       -0.0297      0.0243  
#>          -1  -0.0257     0.0159       -0.0569      0.0055  
#>           0  -0.0148     0.0133       -0.0408      0.0112  
#>           1  -0.0588     0.0210       -0.0999     -0.0177 *
#>           2  -0.1366     0.0381       -0.2112     -0.0619 *
#>           3  -0.1111     0.0344       -0.1786     -0.0437 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## Full QTT curve at selected quantiles
res_qtt <- qdid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res_qtt)
#> 
#> Overall QTT Curve:
#>  Quantile     QTT Std. Error 95% CB Lower 95% CB Upper
#>       0.1  0.0286     0.1840      -0.8675       0.9247
#>       0.2 -0.0545     0.0615      -0.3543       0.2453
#>       0.3 -0.0781     0.0587      -0.3641       0.2079
#>       0.4 -0.0195     0.0442      -0.2346       0.1957
#>       0.5 -0.0614     0.0579      -0.3433       0.2206
#>       0.6 -0.0257     0.0405      -0.2231       0.1717
#>       0.7  0.0046     0.0429      -0.2042       0.2134
#>       0.8 -0.0435     0.0432      -0.2537       0.1668
#>       0.9  0.0000     0.0537      -0.2616       0.2616
#> 
# }
```
