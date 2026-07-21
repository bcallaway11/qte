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
#>  -0.0282        0.0122    -0.0521     -0.0043 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0266     0.0130        0.0012      0.0520 *
#>          -2  -0.0040     0.0122       -0.0279      0.0198  
#>          -1  -0.0247     0.0129       -0.0500      0.0007  
#>           0  -0.0164     0.0112       -0.0385      0.0056  
#>           1  -0.0584     0.0262       -0.1098     -0.0069 *
#>           2  -0.1470     0.0520       -0.2490     -0.0450 *
#>           3  -0.1111     0.0386       -0.1869     -0.0354 *
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
#>       0.1  0.0286     0.1517      -0.6143       0.6715
#>       0.2 -0.0545     0.0573      -0.2974       0.1885
#>       0.3 -0.0912     0.0531      -0.3162       0.1337
#>       0.4 -0.0126     0.0438      -0.1983       0.1731
#>       0.5 -0.0561     0.0485      -0.2618       0.1496
#>       0.6 -0.0257     0.0448      -0.2158       0.1643
#>       0.7  0.0046     0.0423      -0.1746       0.1838
#>       0.8 -0.0435     0.0611      -0.3023       0.2153
#>       0.9  0.0000     0.0490      -0.2079       0.2079
#> 
# }
```
