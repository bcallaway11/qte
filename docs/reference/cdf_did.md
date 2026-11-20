# Difference-in-Differences in Distribution Functions

Computes the Quantile Treatment Effect on the Treated (QTT) and the
Average Treatment Effect on the Treated (ATT) under parallel trends in
the distribution functions of untreated potential outcomes, with
staggered treatment adoption, for panel data or repeated cross sections.

**Identification.** For each group \\g\\ and period \\t\\,
\$\$F\_{Y_t(0)\|G=g}(y) - F\_{Y\_{t-1}(0)\|G=g}(y) = F\_{Y_t(0)\|C}(y) -
F\_{Y\_{t-1}(0)\|C}(y) \quad \text{for all } y,\$\$ with \\C\\ the
comparison group: parallel trends for the distribution function rather
than the mean (DiD) or the quantile function (QDiD). Roth and Sant'Anna
(2023) show this is the condition under which parallel trends holds for
every transformation of the outcome. It point-identifies the
counterfactual distribution, so the QTT, for discrete outcomes too, and
the implied ATT is the usual DiD ATT. Covariates make the assumption
conditional on \\X\\.

See
[`cdf_did_gt`](https://bcallaway11.github.io/qte/reference/cdf_did_gt.md)
for how the counterfactual distribution is estimated (including the
rearrangement into a valid distribution function and the three covariate
methods).

## Usage

``` r
cdf_did(
  yname,
  gname,
  tname,
  idname = NULL,
  data,
  panel = TRUE,
  xformula = ~1,
  est_method = c("aipw", "ipw", "reg"),
  weightsname = NULL,
  clustervars = NULL,
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

- panel:

  Logical. `TRUE` (default) for panel data; `FALSE` for repeated cross
  sections (see the section on repeated cross sections).

- xformula:

  One-sided formula of covariates. Default `~1`. With panel data, each
  unit's covariates are taken from the base period of each group-time
  cell (never from a post-treatment period); with repeated cross
  sections, each observation's own values are used.

- est_method:

  With covariates: `"aipw"` (default; doubly robust), `"ipw"` or `"reg"`
  (outcome regression, as in `did`).

- weightsname:

  Name of the column in `data` containing sampling weights. Default
  `NULL` uses equal weights.

- clustervars:

  Name of a column in `data` to cluster standard errors on, in addition
  to `idname` (as in
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)).
  With the empirical bootstrap, whole clusters are resampled.

- control_group:

  `"notyettreated"` (default) or `"nevertreated"`.

- anticipation:

  Number of time periods before treatment during which units may
  anticipate it. Default 0.

- alp:

  Significance level. Default 0.05.

- cband:

  Logical; uniform confidence bands. Default `TRUE`.

- biters:

  Number of bootstrap iterations. Default 100.

- cl:

  Number of cores for the bootstrap. Default 1.

- gt_type:

  `"att"` (default) for the ATT and event study, or `"qtt"` for the QTT
  curve.

- probs:

  For `gt_type = "qtt"`, the quantiles at which to report the QTT.
  Default `seq(0.05, 0.95, 0.05)`.

## Value

For `gt_type = "att"`, a `pte_emp_boot` object; for `gt_type = "qtt"`, a
`pte_qtt` object. Both include `outcome_discreteness`.

## Discrete outcomes

The counterfactual distribution is identified for discrete outcomes as
well, unlike for rank- or quantile-based estimators such as
[`cic`](https://bcallaway11.github.io/qte/reference/cic.md). Quantile
effects still jump at mass points and their bootstrap inference can be
unreliable, so `gt_type = "qtt"` warns (once per call) when the outcome
looks discrete; see
[`cic`](https://bcallaway11.github.io/qte/reference/cic.md).

## Repeated cross sections

[`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
[`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
[`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md),
`cdf_did()` and `lou_qtt(lagged_outcome_cov = FALSE)` also work with
repeated cross sections (`panel = FALSE`; `idname` not needed). Without
covariates they use only each period's outcome distributions, so on
panel data the estimates are the same either way; the standard errors
differ, because the bootstrap then resamples observations rather than
units. Some things need the same units in several periods and so require
panel data:
[`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) (each
comparison unit's change over time),
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
(the copula between periods), the lagged outcome in
[`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md),
and `gt_type = "qott"` (the distribution of individual treatment
effects).

## References

Roth, Jonathan, and Pedro H. C. Sant'Anna. “When Is Parallel Trends
Sensitive to Functional Form?” Econometrica 91(2), pp. 737-747, 2023.

Sant'Anna, Pedro H. C., and Jun Zhao. “Doubly Robust
Difference-in-Differences Estimators.” Journal of Econometrics 219(1),
pp. 101-122, 2020.

## Examples

``` r
# \donttest{
data(mpdta, package = "did")
res <- cdf_did(yname = "lemp", gname = "first.treat", tname = "year",
               idname = "countyreal", data = mpdta, gt_type = "qtt",
               probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res)
#> 
#> Overall QTT curve:
#>  Quantile     QTT Std. Error [95% Simult. Conf. Band]        
#>       0.1 -0.0403     0.1131                  -0.3415 0.2609 
#>       0.2 -0.0602     0.0596                  -0.2188 0.0985 
#>       0.3 -0.0103     0.0488                  -0.1403 0.1197 
#>       0.4 -0.0690     0.0625                  -0.2354 0.0975 
#>       0.5 -0.0159     0.0436                  -0.1319 0.1002 
#>       0.6 -0.0241     0.0326                  -0.1109 0.0627 
#>       0.7  0.0089     0.0467                  -0.1153 0.1332 
#>       0.8 -0.0192     0.0373                  -0.1184 0.0801 
#>       0.9 -0.0333     0.0438                  -0.1501 0.0834 
#> 
#> Dynamic effects (quantile 0.5):
#>  Event Time     QTT Std. Error [95% Simult. Conf. Band]        
#>          -3  0.0063     0.0792                  -0.1984 0.2111 
#>          -2  0.1155     0.0543                  -0.0384 0.2694 
#>          -1 -0.0609     0.0383                  -0.1833 0.0615 
#>           0 -0.0159     0.0455                  -0.1392 0.1074 
#>           1 -0.0277     0.0526                  -0.1867 0.1313 
#>           2 -0.1320     0.0811                  -0.7681 0.5041 
#>           3 -0.3346     0.0875                  -0.8896 0.2203 
#> 
#> Group effects (quantile 0.5):
#>  Group     QTT Std. Error [95% Simult. Conf. Band]        
#>   2004 -0.1426     0.0612                  -0.7259 0.4408 
#>   2006 -0.0528     0.0422                  -0.2276 0.1221 
#>   2007 -0.1123     0.0850                  -0.4298 0.2053 
#> 
#> ATT, and test of QTT(tau) = ATT for all quantiles:
#>      Parameter     ATT Std. Error [95% Pointwise Conf. Int.] p (QTT = ATT)
#>        Overall -0.0305     0.0097        -0.0494     -0.0115          1.00
#>  Event time -3  0.0298     0.0121         0.0061      0.0534          0.80
#>  Event time -2 -0.0024     0.0107        -0.0234      0.0185          0.35
#>  Event time -1 -0.0243     0.0135        -0.0506      0.0021          0.85
#>   Event time 0 -0.0189     0.0094        -0.0373     -0.0005          1.00
#>   Event time 1 -0.0536     0.0140        -0.0811     -0.0261          0.80
#>   Event time 2 -0.1363     0.0378        -0.2103     -0.0623          0.75
#>   Event time 3 -0.1008     0.0329        -0.1654     -0.0363          0.75
#> 
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Simultaneous bands are uniform over quantiles, for each QTT curve.
#> p (QTT = ATT): bootstrap sup-t test, from the same draws, that QTT(tau) - ATT = 0
#> at every estimated quantile.
# }
```
