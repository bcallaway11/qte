# Panel QTT (Callaway-Li 2019)

Computes Quantile Treatment Effects on the Treated (QTT) and the Average
Treatment Effect on the Treated (ATT) using the three-period panel
identification strategy of Callaway and Li (2019). Handles two-period
(one post + two pre-treatment periods) and staggered treatment adoption
via `ptetools`. Requires panel data for both groups.

**Identification.** Under a copula stability assumption and
distributional parallel trends on changes, the counterfactual
distribution \\F\_{Y(0)^{\text{post}}\|D=1}\\ is recovered from three
periods of panel data without requiring rank invariance.

**pre_copula.** Controls which pre-treatment periods are used as the
copula-transfer base for post-treatment cells (g, g+e):

- `"long"` (default):

  pre2 is as many periods before pre1 as the post period is after it,
  matching the copula window to the event horizon. pre1 is the last
  period at least `anticipation` before g (in units of the time
  variable). Cells where pre2 does not exist in the data are excluded
  from aggregation with two-step weight renormalization.

- `"short"`:

  pre2 is the period just before pre1, always anchored at treatment
  onset.

Pre-test placebo cells always use the sliding base regardless of
`pre_copula`.

## Usage

``` r
panel_qtt(
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
  probs = NULL,
  pre_copula = "long"
)
```

## Arguments

- yname:

  Name of the outcome variable.

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

  One-sided formula for covariates. Default `~1`. Each unit's covariates
  are taken from the earliest period the group-time cell uses (the
  second pre-treatment period), so they are fixed across all three
  periods.

- cond_dist:

  How conditional distributions are estimated when `xformula` has
  covariates: `"qr"` (default) for quantile regression or `"dr"` for
  distribution regression. See
  [`cic_gt`](https://bcallaway11.github.io/qte/reference/cic_gt.md).
  Ignored without covariates.

- weightsname:

  Name of the sampling weights column. Default `NULL`.

- clustervars:

  Name of a column in `data` to cluster standard errors on, in addition
  to `idname` (as in
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)):
  one variable, constant over time within each unit. Default `NULL`
  clusters at the unit level. With the empirical bootstrap, whole
  clusters are resampled.

- control_group:

  `"notyettreated"` (default) or `"nevertreated"`.

- anticipation:

  Anticipation, in units of the time variable. Default `0`.

- alp:

  Significance level. Default `0.05`.

- cband:

  Logical; uniform confidence band if `TRUE` (default).

- biters:

  Number of bootstrap iterations. Default `100`.

- cl:

  Number of parallel clusters. Default `1`.

- gt_type:

  `"att"` (default) or `"qtt"`.

- probs:

  Quantile grid for `gt_type = "qtt"`. Default `seq(0.05, 0.95, 0.05)`.

- pre_copula:

  `"long"` (default) or `"short"`. See Details.

## Value

For `gt_type = "att"`, a `pte_emp_boot` object. For `gt_type = "qtt"`, a
`pte_qtt` object.

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
[`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) (each
comparison unit's change over time), `panel_qtt()` (the copula between
periods), the lagged outcome in
[`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md),
and `gt_type = "qott"` (the distribution of individual treatment
effects).

## References

Callaway, Brantly and Tong Li. “Quantile Treatment Effects in
Difference-in-Differences Models with Panel Data.” Quantitative
Economics 10(4), pp. 1579-1618, 2019.

## See also

[`cic`](https://bcallaway11.github.io/qte/reference/cic.md),
[`qdid`](https://bcallaway11.github.io/qte/reference/qdid.md),
[`mdid`](https://bcallaway11.github.io/qte/reference/mdid.md),
[`ddid`](https://bcallaway11.github.io/qte/reference/ddid.md)

## Examples

``` r
# \donttest{
data(mpdta, package = "did")

## Panel QTT with rolling pre-period copula (default pre_copula = "long")
res <- panel_qtt(yname = "lemp", gname = "first.treat", tname = "year",
                 idname = "countyreal", data = mpdta,
                 gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res)
#> 
#> Overall QTT curve:
#>  Quantile     QTT Std. Error [95% Simult. Conf. Band]        
#>       0.1 -0.0245     0.1423                  -0.4830 0.4341 
#>       0.2 -0.0808     0.0616                  -0.2794 0.1178 
#>       0.3 -0.0011     0.0639                  -0.2072 0.2049 
#>       0.4 -0.0649     0.0618                  -0.2642 0.1344 
#>       0.5 -0.0917     0.0480                  -0.2464 0.0631 
#>       0.6  0.0124     0.0372                  -0.1075 0.1324 
#>       0.7  0.0514     0.0511                  -0.1134 0.2162 
#>       0.8 -0.0065     0.0388                  -0.1317 0.1186 
#>       0.9 -0.0661     0.0317                  -0.1681 0.0360 
#> 
#> Dynamic effects (quantile 0.5):
#>  Event Time     QTT Std. Error [95% Simult. Conf. Band]        
#>          -2  0.0118     0.0771                  -0.2143 0.2378 
#>          -1 -0.0507     0.0531                  -0.2112 0.1098 
#>           0 -0.0442     0.0532                  -0.1926 0.1041 
#>           1 -0.1751     0.0856                  -0.4512 0.1009 
#> 
#> Group effects (quantile 0.5):
#>  Group     QTT Std. Error [95% Simult. Conf. Band]        
#>   2006 -0.1022     0.0621                  -0.2934 0.0890 
#>   2007 -0.1191     0.0585                  -0.3564 0.1182 
#> 
#> ATT, and test of QTT(tau) = ATT for all quantiles:
#>      Parameter     ATT Std. Error [95% Pointwise Conf. Int.] p (QTT = ATT)
#>        Overall -0.0355     0.0120        -0.0590     -0.0120          0.85
#>  Event time -2 -0.0095     0.0218        -0.0521      0.0332          0.75
#>  Event time -1 -0.0415     0.0110        -0.0630     -0.0200          0.30
#>   Event time 0 -0.0301     0.0115        -0.0526     -0.0076          1.00
#>   Event time 1 -0.0612     0.0189        -0.0982     -0.0242          1.00
#> 
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Simultaneous bands are uniform over quantiles, for each QTT curve.
#> p (QTT = ATT): bootstrap sup-t test, from the same draws, that QTT(tau) - ATT = 0
#> at every estimated quantile.
# }
```
