# DiD-Based Quantile Treatment Effects

## Overview

The `qte` package provides six DiD-based estimators for the QTT. All
share the same `yname`/`gname`/`tname`/`idname` interface and support
staggered treatment adoption natively. They differ in their
identification assumptions and the restrictions those assumptions impose
on the data-generating process.

This vignette describes each estimator conceptually and shows the
minimal call. For a full applied workflow with aggregations and
event-study plots, see
[`vignette("staggered-adoption")`](https://bcallaway11.github.io/qte/articles/staggered-adoption.md).

``` r
library(qte)
set.seed(42)
data(mpdta, package = "did")
```

We use `mpdta` throughout: a balanced panel of 500 US counties observed
annually from 2003–2007, with `lemp` (log employment) as the outcome and
`first.treat` as the cohort variable (year of first minimum wage
increase; 0 = never treated).

The common call structure for all six estimators is:

``` r
estimator(
  yname  = "lemp",        # outcome
  gname  = "first.treat", # first treatment period (0 = never treated)
  tname  = "year",        # calendar time
  idname = "countyreal",  # unit id
  data   = mpdta,
  gt_type = "att",        # "att" or "qtt"
  biters  = 100
)
```

------------------------------------------------------------------------

## Change in Changes — `cic()`

**Reference:** Athey and Imbens ([2006](#ref-athey-imbens-2006))

**Identification assumption:** The rank of each unit in the untreated
outcome distribution is stable over time (rank invariance), and the
support of the treated group’s pre-period distribution is contained in
the support of the untreated group’s pre-period distribution (support
condition). No functional form for the outcome equation is assumed.

**When to use:** A flexible nonparametric alternative to linear DiD when
the parallel trends assumption is implausible in levels but rank
stability is credible. Works with panel or repeated cross sections.

**Key restriction:** Stronger than standard parallel trends — requires
the entire rank ordering to be preserved, not just the mean.

``` r
res_cic <- cic(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 50
)
summary(res_cic)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.] 
#>  -0.0197        0.0185    -0.0532      0.0139 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0508     0.0176        0.0164      0.0853 *
#>          -2   0.0158     0.0159       -0.0154      0.0471  
#>          -1  -0.0128     0.0156       -0.0434      0.0178  
#>           0  -0.0081     0.0185       -0.0442      0.0281  
#>           1  -0.0364     0.0196       -0.0748      0.0020  
#>           2  -0.1226     0.0348       -0.1907     -0.0544 *
#>           3  -0.0930     0.0369       -0.1653     -0.0206 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

------------------------------------------------------------------------

## Quantile Difference-in-Differences — `qdid()`

**Reference:** Callaway et al. ([2018](#ref-callaway-li-oka-2018))

**Identification assumption:** Distributional parallel trends — the
entire untreated potential outcome distribution shifts in parallel
across time for treated and untreated groups (not just the mean).
Formally, $F_{Y(0)|D = 1,t}$ shifts over time in the same way as
$F_{Y(0)|D = 0,t}$.

**When to use:** A direct distributional analog of standard DiD. Weaker
than CiC in that it does not require rank invariance, but stronger than
mean DiD. Works with panel or repeated cross sections.

**Key restriction:** Requires the full distributional shift to be
parallel, not just the mean shift; this is testable in pre-treatment
periods.

``` r
res_qdid <- qdid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 50
)
summary(res_qdid)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.] 
#>  -0.0271        0.0136    -0.0567      0.0026 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0266     0.0145       -0.0018      0.0550  
#>          -2  -0.0027     0.0143       -0.0308      0.0254  
#>          -1  -0.0257     0.0143       -0.0538      0.0023  
#>           0  -0.0148     0.0129       -0.0401      0.0105  
#>           1  -0.0588     0.0196       -0.0972     -0.0203 *
#>           2  -0.1366     0.0374       -0.2099     -0.0632 *
#>           3  -0.1111     0.0400       -0.1895     -0.0328 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

------------------------------------------------------------------------

## Distributional Difference-in-Differences — `ddid()`

**Reference:** Callaway et al. ([2018](#ref-callaway-li-oka-2018))

**Identification assumption:** Conditional distributional parallel
trends — the untreated potential outcome distribution shifts in parallel
conditional on pre-treatment covariates. Extends QDiD to allow covariate
adjustment.

**When to use:** When distributional parallel trends holds only after
conditioning on covariates. The covariate adjustment is done via outcome
regression on the CDFs, making it more flexible than linear regression
adjustment for standard DiD.

**Key restriction:** Same spirit as QDiD but allows covariates. Requires
panel data.

``` r
res_ddid <- ddid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 50
)
summary(res_ddid)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -0.0453        0.0139    -0.0719     -0.0187 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0134     0.0144       -0.0148      0.0416  
#>          -2  -0.0146     0.0113       -0.0367      0.0075  
#>          -1  -0.0403     0.0157       -0.0711     -0.0096 *
#>           0  -0.0324     0.0128       -0.0575     -0.0072 *
#>           1  -0.0637     0.0204       -0.1038     -0.0237 *
#>           2  -0.1377     0.0415       -0.2192     -0.0563 *
#>           3  -0.1087     0.0393       -0.1856     -0.0317 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

------------------------------------------------------------------------

## Mean Difference-in-Differences — `mdid()`

**Reference:** Callaway and Li ([2019](#ref-callaway-li-2019))

**Identification assumption:** Standard (mean) parallel trends — the
mean of the untreated potential outcome evolves in parallel for treated
and untreated groups. This is the same assumption as in standard DiD.

**When to use:** When only mean parallel trends is credible, not
distributional parallel trends.
[`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) recovers
the ATT under this weaker assumption and can also return the QTT by
combining mean DiD with a rank-based counterfactual distribution.

**Key restriction:** Weakest distributional assumption of the four DiD
methods. The QTT estimate relies on additional rank-stability within the
mean DiD framework.

``` r
res_mdid <- mdid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 50
)
summary(res_mdid)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -0.0305        0.0109    -0.0565     -0.0044 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0298     0.0148        0.0007      0.0588 *
#>          -2  -0.0024     0.0116       -0.0252      0.0203  
#>          -1  -0.0243     0.0150       -0.0537      0.0051  
#>           0  -0.0189     0.0100       -0.0386      0.0007  
#>           1  -0.0536     0.0159       -0.0847     -0.0225 *
#>           2  -0.1363     0.0373       -0.2093     -0.0632 *
#>           3  -0.1008     0.0315       -0.1625     -0.0391 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

------------------------------------------------------------------------

## Panel QTT — `panel_qtt()`

**Reference:** Callaway and Li ([2019](#ref-callaway-li-2019))

**Identification assumption:** Copula stability — the dependence
structure (copula) between the pre-treatment and post-treatment
untreated potential outcomes is the same for treated and untreated
groups. This is a panel-data assumption that uses three periods: two
pre-treatment periods to estimate the copula, and one post-treatment
period.

**When to use:** When panel data with at least two pre-treatment periods
per cohort are available and the copula stability assumption is
credible. Uniquely identifies the individual-level rank mapping from pre
to post, enabling inference on the QTT under weaker distributional
assumptions than CiC.

**Key restriction:** Requires panel data with sufficient pre-treatment
periods. The `pre_copula` argument controls whether a fixed (`"short"`)
or rolling (`"long"`, default) base period is used.

``` r
res_pqtt <- panel_qtt(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 50
)
summary(res_pqtt)
#> 
#> Overall QTT Curve:
#>  Quantile     QTT Std. Error 95% CB Lower 95% CB Upper
#>       0.1 -0.0270     0.1296      -0.5674       0.5133
#>       0.2 -0.0852     0.0711      -0.3818       0.2115
#>       0.3 -0.0250     0.0601      -0.2756       0.2257
#>       0.4 -0.0680     0.0553      -0.2984       0.1625
#>       0.5 -0.0965     0.0530      -0.3176       0.1246
#>       0.6  0.0070     0.0385      -0.1536       0.1675
#>       0.7  0.0409     0.0622      -0.2187       0.3004
#>       0.8 -0.0201     0.0836      -0.3686       0.3284
#>       0.9 -0.0832     0.0603      -0.3344       0.1680
```

------------------------------------------------------------------------

## Lagged-Outcome Unconfoundedness — `lou_qte()`

**Identification assumption:** Unconfoundedness conditional on the
lagged outcome — treatment is as good as randomly assigned after
conditioning on the pre-treatment value of the outcome. This is a
selection-on-observables assumption, not a parallel trends assumption.

**When to use:** When the pre-treatment outcome is a strong proxy for
the unobserved confounders driving selection into treatment, and the
parallel trends assumption is implausible. Common in labor and health
economics where lagged earnings or health status are strong predictors
of treatment.

**Key restriction:** Does not require parallel trends, but assumes that
conditioning on the lagged outcome absorbs all confounding. Requires
panel data.

``` r
res_lou <- lou_qte(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 50
)
summary(res_lou)
#> 
#> Overall ATT:  
#>      ATT    Std. Error     [ 95%  Conf. Int.] 
#>  -0.0384        0.0221    -0.0831      0.0063 
#> 
#> 
#> Dynamic Effects:
#>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -3   0.0294     0.0135        0.0031      0.0558 *
#>          -2  -0.0124     0.0149       -0.0417      0.0168  
#>          -1  -0.0310     0.0189       -0.0680      0.0060  
#>           0  -0.0243     0.0206       -0.0646      0.0160  
#>           1  -0.0869     0.0356       -0.1568     -0.0171 *
#>           2  -0.1394     0.0328       -0.2036     -0.0752 *
#>           3  -0.1093     0.0386       -0.1850     -0.0336 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

------------------------------------------------------------------------

## Summary

| Function                                                                  | Assumption                      | Panel required | Notes                               |
|---------------------------------------------------------------------------|---------------------------------|----------------|-------------------------------------|
| [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md)             | Rank invariance + support       | Optional       | Nonparametric; strongest assumption |
| [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md)           | Distributional parallel trends  | Optional       | DiD analog for distributions        |
| [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md)           | Conditional distributional PT   | Yes            | Adds covariate adjustment to QDiD   |
| [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md)           | Mean parallel trends            | Optional       | Weakest DiD assumption              |
| [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md) | Copula stability                | Yes            | Requires ≥ 2 pre-periods per cohort |
| [`lou_qte()`](https://bcallaway11.github.io/qte/reference/lou_qte.md)     | Lagged-outcome unconfoundedness | Yes            | Not a DiD assumption                |

All estimators return either a `pte_results` object (`gt_type = "att"`)
or a `pte_qtt` object (`gt_type = "qtt"`), and support `autoplot()` for
plotting overall curves and event-study summaries. See
[`vignette("staggered-adoption")`](https://bcallaway11.github.io/qte/articles/staggered-adoption.md)
for the full applied workflow.

## References

Athey, Susan, and Guido Imbens. 2006. “Identification and Inference in
Nonlinear Difference-in-Differences Models.” *Econometrica* 74 (2):
431–97.

Callaway, Brantly, and Tong Li. 2019. “Quantile Treatment Effects in
Difference in Differences Models with Panel Data.” *Quantitative
Economics* 10 (4): 1579–618.

Callaway, Brantly, Tong Li, and Tatsushi Oka. 2018. “Quantile Treatment
Effects in Difference in Differences Models Under Dependence
Restrictions and with Only Two Time Periods.” *Journal of Econometrics*
206 (2): 395–413.
