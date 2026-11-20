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

  One-sided formula for covariates. Default `~1` uses no covariates.
  With panel data, each unit's covariates are taken from the base period
  of each group-time cell (never from a post-treatment period); with
  repeated cross sections, each observation's own values are used.

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

[`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
[`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
`mdid()`,
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
#>  -0.0305        0.0128    -0.0559      -0.005 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0298     0.0127        0.0049      0.0546 *
#>          -2  -0.0024     0.0122       -0.0263      0.0214  
#>          -1  -0.0243     0.0129       -0.0495      0.0009  
#>           0  -0.0189     0.0122       -0.0428      0.0050  
#>           1  -0.0536     0.0110       -0.0752     -0.0320 *
#>           2  -0.1363     0.0449       -0.2243     -0.0482 *
#>           3  -0.1008     0.0359       -0.1712     -0.0304 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 

## Full QTT curve at selected quantiles
res_qtt <- mdid(yname = "lemp", gname = "first.treat", tname = "year",
                idname = "countyreal", data = mpdta,
                gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res_qtt)
#> 
#> Overall QTT curve:
#>  Quantile     QTT Std. Error [95% Simult. Conf. Band]        
#>       0.1 -0.0125     0.1055                  -0.2604 0.2354 
#>       0.2 -0.0452     0.0888                  -0.2539 0.1634 
#>       0.3 -0.0118     0.0466                  -0.1213 0.0977 
#>       0.4 -0.0760     0.0404                  -0.1710 0.0191 
#>       0.5 -0.0429     0.0493                  -0.1588 0.0730 
#>       0.6 -0.0260     0.0294                  -0.0952 0.0432 
#>       0.7  0.0088     0.0413                  -0.0883 0.1058 
#>       0.8 -0.0071     0.0290                  -0.0753 0.0611 
#>       0.9 -0.0144     0.0455                  -0.1213 0.0925 
#> 
#> Dynamic effects (quantile 0.5):
#>  Event Time     QTT Std. Error [95% Simult. Conf. Band]        
#>          -3  0.0338     0.0311                  -0.0705 0.1381 
#>          -2  0.0510     0.0431                  -0.0913 0.1932 
#>          -1 -0.0198     0.0635                  -0.2200 0.1803 
#>           0 -0.0407     0.0468                  -0.1620 0.0805 
#>           1 -0.0685     0.0561                  -0.2447 0.1077 
#>           2 -0.0866     0.1125                  -0.4029 0.2297 
#>           3 -0.1680     0.0944                  -0.5182 0.1823 
#> 
#> Group effects (quantile 0.5):
#>  Group     QTT Std. Error [95% Simult. Conf. Band]        
#>   2004 -0.1098     0.0748                  -0.3268 0.1071 
#>   2006 -0.0178     0.0590                  -0.1865 0.1508 
#>   2007 -0.1046     0.0759                  -0.3624 0.1533 
#> 
#> ATT, and test of QTT(tau) = ATT for all quantiles:
#>      Parameter     ATT Std. Error [95% Pointwise Conf. Int.] p (QTT = ATT)
#>        Overall -0.0305     0.0116        -0.0531     -0.0078          0.95
#>  Event time -3  0.0298     0.0125         0.0054      0.0542          0.85
#>  Event time -2 -0.0024     0.0138        -0.0294      0.0245          0.85
#>  Event time -1 -0.0243     0.0136        -0.0509      0.0023          0.95
#>   Event time 0 -0.0189     0.0126        -0.0435      0.0057          1.00
#>   Event time 1 -0.0536     0.0172        -0.0872     -0.0199          0.35
#>   Event time 2 -0.1363     0.0388        -0.2123     -0.0602          0.10
#>   Event time 3 -0.1008     0.0361        -0.1717     -0.0300          0.60
#> 
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Simultaneous bands are uniform over quantiles, for each QTT curve.
#> p (QTT = ATT): bootstrap sup-t test, from the same draws, that QTT(tau) - ATT = 0
#> at every estimated quantile.
# }
```
