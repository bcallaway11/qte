# Staggered Treatment Adoption with the qte Package

## Introduction

In many empirical settings, units adopt treatment at different points in
time — a structure called *staggered treatment adoption*. For example,
states may adopt a policy in different years, or firms may receive an
intervention at different times during a study period.

All seven panel estimators in `qte` support staggered adoption via the
same interface: specify the outcome (`yname`), the
first-treatment-period cohort variable (`gname`, with 0 for
never-treated units), calendar time (`tname`), and a unit identifier
(`idname`). The package follows the approach of Callaway and Sant’Anna
([2021](#ref-callaway-santanna-2021)), computing group-time effects
$`QTT(g, t)`$ for every cohort $`g`$ and period $`t`$, then aggregating
these into overall, group-specific, and event-study (dynamic) summaries.

The primary output of interest is the **QTT curve** — the treatment
effect at each quantile $`\tau`$ of the treated group’s outcome
distribution. Each QTT run also reports the ATT implied by the same
assumption, and tests whether the QTT curve is flat at the ATT.

This vignette demonstrates the full applied workflow using the `mpdta`
dataset. For conceptual background on the identification assumptions
behind each estimator, see
[`vignette("panel-estimators")`](https://bcallaway11.github.io/qte/articles/panel-estimators.md).

``` r
library(qte)
library(ggplot2)
library(ptetools)
set.seed(42)
data(mpdta, package = "did")
```

`mpdta` is a balanced panel of 500 US counties observed annually from
2003 to 2007. The outcome is `lemp` (log county-level employment).
Counties are grouped by `first.treat`, the year they first adopted a
minimum wage increase (2004, 2006, or 2007); counties with
`first.treat == 0` never adopt.

``` r
table(mpdta$first.treat)
#> 
#>    0 2004 2006 2007 
#> 1545  100  200  655
```

------------------------------------------------------------------------

## QTT curve estimation

We begin with
[`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) using
`gt_type = "qtt"` to estimate the full QTT curve. For each aggregation,
the group-time cells’ distribution functions are mixed first (with the
same weights as the corresponding ATT aggregation) and the mixture is
then inverted at each requested quantile. This gives the QTT for the
treated population as a whole, which averaging cell-level QTTs would
not.

``` r
res_qtt <- cic(
  yname   = "lemp",
  gname   = "first.treat",
  tname   = "year",
  idname  = "countyreal",
  data    = mpdta,
  biters  = 100,
  gt_type = "qtt",
  probs   = seq(0.1, 0.9, by = 0.1)
)
```

[`summary()`](https://rdrr.io/r/base/summary.html) prints the overall
QTT curve, the event-study and cohort QTTs at the median (change this
with `probs`), and the ATTs:

``` r
summary(res_qtt)
#> 
#> Overall QTT curve:
#>  Quantile     QTT Std. Error [95% Simult. Conf. Band]        
#>       0.1  0.0278     0.0932                  -0.2339 0.2895 
#>       0.2 -0.0230     0.0749                  -0.2331 0.1871 
#>       0.3 -0.0103     0.0466                  -0.1412 0.1206 
#>       0.4 -0.0690     0.0534                  -0.2188 0.0808 
#>       0.5 -0.0185     0.0558                  -0.1752 0.1382 
#>       0.6 -0.0149     0.0412                  -0.1305 0.1007 
#>       0.7  0.0389     0.0602                  -0.1301 0.2078 
#>       0.8 -0.0387     0.0505                  -0.1803 0.1030 
#>       0.9  0.0135     0.0934                  -0.2486 0.2755 
#> 
#> Dynamic effects (quantile 0.5):
#>  Event Time     QTT Std. Error [95% Simult. Conf. Band]        
#>          -3  0.0159     0.0653                  -0.2166 0.2485 
#>          -2  0.0868     0.0792                  -0.2502 0.4238 
#>          -1 -0.0683     0.0601                  -0.2492 0.1126 
#>           0 -0.0237     0.0572                  -0.2319 0.1844 
#>           1 -0.0429     0.0609                  -0.2154 0.1295 
#>           2 -0.0960     0.1034                  -0.4330 0.2409 
#>           3 -0.1972     0.1190                  -0.6089 0.2145 
#> 
#> Group effects (quantile 0.5):
#>  Group     QTT Std. Error [95% Simult. Conf. Band]        
#>   2004 -0.1391     0.0908                  -0.4426 0.1645 
#>   2006  0.0149     0.0591                  -0.1809 0.2107 
#>   2007 -0.1181     0.0691                  -0.3238 0.0875 
#> 
#> ATT, and test of QTT(tau) = ATT for all quantiles:
#>      Parameter     ATT Std. Error [95% Pointwise Conf. Int.] p (QTT = ATT)
#>        Overall -0.0197     0.0143        -0.0477      0.0083          0.96
#>  Event time -3  0.0508     0.0244         0.0030      0.0987          0.97
#>  Event time -2  0.0158     0.0138        -0.0113      0.0430          0.76
#>  Event time -1 -0.0128     0.0138        -0.0399      0.0142          0.85
#>   Event time 0 -0.0081     0.0140        -0.0354      0.0193          0.97
#>   Event time 1 -0.0364     0.0228        -0.0811      0.0083          0.91
#>   Event time 2 -0.1226     0.0405        -0.2019     -0.0433          0.34
#>   Event time 3 -0.0930     0.0413        -0.1739     -0.0120          0.42
#> 
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> Simultaneous bands are uniform over quantiles, for each QTT curve.
#> p (QTT = ATT): bootstrap sup-t test, from the same draws, that QTT(tau) - ATT = 0
#> at every estimated quantile.
```

### Overall QTT curve

[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
plots the overall curve with a uniform confidence band, and the overall
ATT as a horizontal line:

``` r
autoplot(res_qtt)
```

![Overall QTT curve from cic(), with the overall
ATT](staggered-adoption_files/figure-html/qtt-plot-1.png)

### Dynamic QTT — event study by quantile

`type = "dynamic"` plots the QTT at selected quantile(s) across event
times. The `plot_probs` argument selects which quantiles to show.

``` r
autoplot(res_qtt, type = "dynamic", plot_probs = 0.5)
```

![Event-study plot for the median
QTT](staggered-adoption_files/figure-html/qtt-dynamic-single-1.png)

Multiple quantiles can be overlaid by passing a vector to `plot_probs`
(values must be in the estimated `probs` grid). Confidence intervals are
then hidden by default (`plot_ci = TRUE` shows them):

``` r
autoplot(res_qtt, type = "dynamic", plot_probs = c(0.1, 0.5, 0.9))
```

![Event-study plot for multiple quantiles of the
QTT](staggered-adoption_files/figure-html/qtt-dynamic-multi-1.png)

Pre-treatment estimates near zero support the identifying assumption.
The post-treatment pattern shows the estimated QTT in each period after
adoption.

### QTT by cohort

`type = "group"` shows each cohort’s QTT, averaged over its
post-treatment periods:

``` r
autoplot(res_qtt, type = "group", plot_probs = 0.5)
```

![Median QTT by treatment
cohort](staggered-adoption_files/figure-html/qtt-group-1.png)

------------------------------------------------------------------------

## Is the effect more than a location shift?

If treatment shifts the whole distribution by the same amount, the QTT
equals the ATT at every quantile. Because the QTT run also computes the
ATT from the same bootstrap draws, the difference $`QTT(\tau) - ATT`$
has its own uniform confidence band, and the last column of the ATT
table in [`summary()`](https://rdrr.io/r/base/summary.html) above is the
p-value of the test that $`QTT(\tau) = ATT`$ at every estimated
quantile. This is a valid test, unlike checking whether the ATT lies
inside the QTT band, which ignores that the two estimates are
correlated.

``` r
autoplot(res_qtt, type = "qtt_minus_att")
```

![QTT minus ATT across quantiles, with a uniform confidence
band](staggered-adoption_files/figure-html/qtt-minus-att-1.png)

Here the band covers zero at every quantile (p-value 0.96), so there is
no evidence that the minimum wage changed the shape of the employment
distribution beyond shifting it. In general, a band that excludes zero
at some quantile rejects a pure location shift: the treatment changes
the shape of the outcome distribution, which the ATT alone would miss.
Note that this compares distributions. Even a location shift of the
distribution is consistent with individual effects that vary across
units.

------------------------------------------------------------------------

## ATT and event-study summary

When only the ATT is needed, `gt_type = "att"` is faster: it skips the
distribution functions. Event times $`e < 0`$ are pre-treatment periods;
$`e \geq 0`$ are post-treatment.

``` r
res_att <- cic(
  yname   = "lemp",
  gname   = "first.treat",
  tname   = "year",
  idname  = "countyreal",
  data    = mpdta,
  biters  = 100,
  gt_type = "att"
)
summary(res_att)
```

    #>
    #> Overall ATT:
    #>      ATT    Std. Error     [ 95%  Conf. Int.]
    #>  -0.0197        0.0149     -0.052      0.0127
    #>
    #>
    #> Dynamic Effects:
    #>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]
    #>          -3   0.0508     0.0230        0.0057      0.0960 *
    #>          -2   0.0158     0.0151       -0.0137      0.0454
    #>          -1  -0.0128     0.0146       -0.0414      0.0158
    #>           0  -0.0081     0.0136       -0.0348      0.0187
    #>           1  -0.0364     0.0251       -0.0855      0.0127
    #>           2  -0.1226     0.0445       -0.2098     -0.0353 *
    #>           3  -0.0930     0.0495       -0.1901      0.0041
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

``` r
autoplot(res_att, type = "dynamic")
```

![Event-study plot of ATT estimates by event
time](staggered-adoption_files/figure-html/att-event-study-1.png)

------------------------------------------------------------------------

## Comparing estimators

All seven estimators share the same interface. Here we run each with
`gt_type = "qtt"` and `probs = 0.5` to compare median QTT estimates. The
estimators are listed in the same order as in
[`vignette("panel-estimators")`](https://bcallaway11.github.io/qte/articles/panel-estimators.md).

``` r
res_cic <- cic(
  yname = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = 0.5, biters = 100
)
# ... and the same call for qdid(), panel_qtt(), ddid(), mdid(), cdf_did()
# and lou_qtt(), saved as res_qdid, res_pqtt, res_ddid, res_mdid, res_cdf
# and res_lou
```

``` r
median_qtt <- function(res, label) {
  row <- res$overall[res$overall$probs == 0.5, ]
  data.frame(estimator = label, qtt = round(row$qtt, 4), se = round(row$se, 4),
             att = round(res$att$overall$att, 4))
}

do.call(rbind, list(
  median_qtt(res_cic,  "cic()"),
  median_qtt(res_qdid, "qdid()"),
  median_qtt(res_pqtt, "panel_qtt()"),
  median_qtt(res_ddid, "ddid()"),
  median_qtt(res_mdid, "mdid()"),
  median_qtt(res_cdf,  "cdf_did()"),
  median_qtt(res_lou,  "lou_qtt()")
))
#>     estimator     qtt     se     att
#> 1       cic() -0.0185 0.0557 -0.0197
#> 2      qdid() -0.0406 0.0483 -0.0282
#> 3 panel_qtt() -0.0917 0.0549 -0.0358
#> 4      ddid() -0.0884 0.0487 -0.0453
#> 5      mdid() -0.0429 0.0442 -0.0305
#> 6   cdf_did() -0.0159 0.0556 -0.0305
#> 7   lou_qtt()  0.0107 0.1041 -0.0384
```

Estimates will generally differ because each estimator relies on a
different identifying assumption. Similarity across methods strengthens
confidence in the findings; divergence warrants investigation of which
assumption is more plausible for the application. The ATTs differ too,
since each assumption implies its own counterfactual mean.
[`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) and
[`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
both imply mean parallel trends, so their ATT is the usual DiD estimate.
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)’s
default long copula also drops group-time cells without enough
pre-treatment periods, so it averages over fewer cells.

------------------------------------------------------------------------

## Covariates

Each estimator’s assumption can be made conditional on covariates
through `xformula`, as with
[`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html).
How they enter depends on the estimator:

- [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  and [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md)
  estimate conditional outcome distributions, by quantile regression
  (`cond_dist = "qr"`, the default) or distribution regression
  (`cond_dist = "dr"`, a logit of $`1\{Y \le y\}`$ on the covariates at
  a grid of thresholds). Distribution regression also handles outcomes
  with mass points.
- [`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
  and
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  take `est_method = c("aipw", "ipw", "reg")`: doubly robust, inverse
  probability weighting, or outcome regression
  ([`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)’s
  outcome regression also takes `cond_dist`).
- [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md)
  adjusts the comparison group’s mean change by a linear regression on
  the covariates.

With panel data, covariates that change over time are taken from each
group-time cell’s base period (for
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md),
the earliest period the cell uses), so the estimators never condition on
covariates measured after treatment. If treated units in some cell have
covariate values the comparison group doesn’t span (say, a region with
treated units only), the estimators warn, listing the cells;
[`ptetools::check_overlap()`](https://rdrr.io/pkg/ptetools/man/check_overlap.html)
reports overlap diagnostics for every cell.

Here we condition on log county population, `lpop`:

``` r
res_cdf_x <- cdf_did(
  yname      = "lemp",
  gname      = "first.treat",
  tname      = "year",
  idname     = "countyreal",
  data       = mpdta,
  xformula   = ~lpop,
  est_method = "aipw",
  biters     = 100,
  gt_type    = "qtt",
  probs      = seq(0.1, 0.9, by = 0.1)
)

res_cic_x <- cic(
  yname     = "lemp",
  gname     = "first.treat",
  tname     = "year",
  idname    = "countyreal",
  data      = mpdta,
  xformula  = ~lpop,
  cond_dist = "dr",
  biters    = 100,
  gt_type   = "qtt",
  probs     = seq(0.1, 0.9, by = 0.1)
)
```

``` r
autoplot(res_cdf_x)
```

![Overall QTT curve from cdf_did() with
covariates](staggered-adoption_files/figure-html/covariates-plot-cdf-1.png)

``` r
autoplot(res_cic_x)
```

![Overall QTT curve from cic() with
covariates](staggered-adoption_files/figure-html/covariates-plot-cic-1.png)

------------------------------------------------------------------------

## Repeated cross sections

If the data are repeated cross sections rather than a panel, set
`panel = FALSE` and omit `idname`.
[`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
[`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
[`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) and
[`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
support repeated cross sections, as does
[`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
with `lagged_outcome_cov = FALSE`.
[`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) and
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
need panel data, as does the distribution of individual effects
(`gt_type = "qott"`).

``` r
res_rcs <- cic(
  yname   = "lemp",
  gname   = "first.treat",
  tname   = "year",
  data    = mpdta,   # no idname
  panel   = FALSE,
  gt_type = "qtt",
  biters  = 100
)
```

------------------------------------------------------------------------

## Inference notes

Standard errors are computed with the empirical bootstrap, resampling
units (panel) or observations (repeated cross sections) with
replacement. The `biters` argument controls the number of iterations
(default 100). Parallel computation is available via the `cl` argument
(`cl = 4` uses 4 cores).

- **Clustering.** `clustervars` resamples whole clusters at a level
  above the unit (for example, states), as in
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html).
- **Sampling weights.** `weightsname` names a column of sampling
  weights, used in every estimation step and in the aggregation across
  cohorts.
- **Confidence bands.** With `cband = TRUE` (the default) each QTT curve
  and event study gets a uniform band, with simultaneous coverage over
  all quantiles or event times. Set `cband = FALSE` for pointwise
  intervals;
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  also takes `cband` to choose which to draw.
- **Discrete outcomes.** The estimators warn when the outcome looks
  discrete (a large mass point, or few distinct values). Quantile
  treatment effects are then usually not point identified, and bootstrap
  inference can be unreliable;
  [`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
  is the exception for identification.

## References

Callaway, Brantly, and Pedro HC Sant’Anna. 2021.
“Difference-in-Differences with Multiple Time Periods.” *Journal of
Econometrics* 225 (2): 200–230.
