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

- `"or"`:

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
lou_qte(
  yname,
  gname,
  tname,
  idname = NULL,
  data,
  xformula = ~1,
  lagged_outcome_cov = TRUE,
  est_method = c("ipw", "or", "aipw"),
  panel = TRUE,
  weightsname = NULL,
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
  beyond the lagged outcome (when `lagged_outcome_cov = TRUE`).

- lagged_outcome_cov:

  Logical; if `TRUE` (default), the pre-treatment outcome is added as a
  covariate in both the propensity score and outcome models. Requires
  `panel = TRUE`.

- est_method:

  Estimation method: `"ipw"` (default), `"or"`, or `"aipw"`. See
  Details.

- panel:

  Logical; `TRUE` (default) for panel data, `FALSE` for repeated cross
  sections. `lagged_outcome_cov = TRUE` requires `panel = TRUE`.

- weightsname:

  Name of the column in `data` containing sampling weights. Default
  `NULL` uses equal weights.

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
res_att <- lou_qte(
  yname = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  lagged_outcome_cov = TRUE, est_method = "ipw",
  gt_type = "att", biters = 20
)
summary(res_att)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -0.0384        0.0132    -0.0599     -0.0169 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0294     0.0121        0.0056      0.0532 *
#>          -2  -0.0124     0.0100       -0.0321      0.0072  
#>          -1  -0.0310     0.0172       -0.0647      0.0028  
#>           0  -0.0243     0.0125       -0.0487      0.0001  
#>           1  -0.0869     0.0305       -0.1467     -0.0271 *
#>           2  -0.1394     0.0314       -0.2010     -0.0778 *
#>           3  -0.1093     0.0400       -0.1877     -0.0309 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## QTT with doubly-robust estimation
res_qtt <- lou_qte(
  yname = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  lagged_outcome_cov = TRUE, est_method = "aipw",
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20
)
summary(res_qtt)
#> 
#> Overall QTT Curve:
#>  Quantile     QTT Std. Error 95% CB Lower 95% CB Upper
#>       0.1  0.0883     0.0747      -0.2475       0.4241
#>       0.2 -0.0315     0.0700      -0.3462       0.2832
#>       0.3 -0.0465     0.0577      -0.3058       0.2128
#>       0.4 -0.0291     0.0507      -0.2569       0.1987
#>       0.5 -0.0330     0.0465      -0.2419       0.1759
#>       0.6 -0.0362     0.0378      -0.2062       0.1338
#>       0.7 -0.0089     0.0437      -0.2053       0.1876
#>       0.8 -0.0187     0.0843      -0.3976       0.3601
#>       0.9 -0.0372     0.0906      -0.4444       0.3700
#> 
# }
```
