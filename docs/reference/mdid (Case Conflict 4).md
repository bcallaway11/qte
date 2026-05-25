# Mean Difference-in-Differences

Computes Quantile Treatment effects on the Treated (QTT) and the Average
Treatment Effect on the Treated (ATT) using the Mean
Difference-in-Differences identification strategy. Handles two-period
data and staggered treatment adoption uniformly. Supports both panel and
repeated cross sections data.

**Identification.** MDiD assumes the counterfactual distribution of
untreated potential outcomes is a location shift of the treated group's
pre-treatment distribution. The size of the shift is the mean DiD
\\\Delta = E\[Y\_{\text{post}}\|D=0\] - E\[Y\_{\text{pre}}\|D=0\]\\, so
\\Q\_{Y(0),\text{post}\|D=1}(\tau) = Q\_{Y,\text{pre}\|D=1}(\tau) +
\Delta\\. This is stronger than parallel trends in means alone: it
additionally requires that the shape of the treated group's outcome
distribution is unchanged in the counterfactual. MDiD is a special case
of QDiD that applies a single mean shift rather than a rank-specific
distributional shift.

**Covariate adjustment.** When `xformula` is specified, the scalar shift
is replaced by a unit-specific conditional mean shift \\\Delta(X_i) =
E\[Y\_{\text{post}}\|D=0,X_i\] - E\[Y\_{\text{pre}}\|D=0,X_i\]\\,
estimated by weighted OLS. The counterfactual for treated unit \\i\\ is
\\Y\_{\text{pre},i} + \Delta(X_i)\\, and the unconditional
counterfactual distribution is the empirical CDF of these values. By the
law of iterated expectations, this consistently estimates
\\F\_{Y(0),\text{post}\|D=1}\\.

## Usage

``` r
mdid(
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

Thuysbaert, Bram. “Distributional Comparisons in Difference in
Differences Models.” Working Paper, 2007.

## Examples

``` r
# \donttest{
data(mpdta, package = "did")

## ATT aggregated across all groups and periods
res_att <- mdid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "att", biters = 20)
summary(res_att)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -0.0305         0.013    -0.0535     -0.0074 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0298     0.0140        0.0024      0.0571 *
#>          -2  -0.0024     0.0139       -0.0296      0.0247  
#>          -1  -0.0243     0.0170       -0.0576      0.0090  
#>           0  -0.0189     0.0124       -0.0433      0.0055  
#>           1  -0.0536     0.0140       -0.0811     -0.0261 *
#>           2  -0.1363     0.0317       -0.1984     -0.0742 *
#>           3  -0.1008     0.0295       -0.1587     -0.0430 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## Full QTT curve at selected quantiles
res_qtt <- mdid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res_qtt)
#> 
#> Overall QTT Curve:
#>  Quantile     QTT Std. Error 95% CB Lower 95% CB Upper
#>       0.1  0.0474     0.1146      -0.4434       0.5381
#>       0.2 -0.0724     0.0783      -0.4077       0.2629
#>       0.3 -0.0465     0.0251      -0.1539       0.0609
#>       0.4 -0.0195     0.0373      -0.1792       0.1402
#>       0.5 -0.0384     0.0415      -0.2162       0.1395
#>       0.6 -0.0305     0.0418      -0.2093       0.1483
#>       0.7 -0.0301     0.0353      -0.1812       0.1211
#>       0.8 -0.0540     0.0422      -0.2348       0.1268
#>       0.9 -0.0185     0.0387      -0.1841       0.1471
#> 
# }
```
