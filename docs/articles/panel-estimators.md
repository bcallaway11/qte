# Panel Data Estimators for Quantile Treatment Effects

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
library(ptetools)
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

**Identifying assumption:** The distributional time-change is the same
for both groups:
$$Q_{Y_{t}(0)|D = 1}\!\left( F_{Y_{t - 1}(0)|D = 1}(y) \right) = Q_{Y_{t}(0)|D = 0}\!\left( F_{Y_{t - 1}(0)|D = 0}(y) \right)\quad{\text{for all}\mspace{6mu}}y,$$
along with a support condition:
$\text{supp}(Y_{t - 1}(0)|D = 1) \subseteq \text{supp}(Y_{t - 1}(0)|D = 0)$.

**Intuition:** The expression
$Q_{Y_{t}(0)|D = g}(F_{Y_{t - 1}(0)|D = g}(y))$ takes a pre-period value
$y$, finds its quantile rank in group $g$’s pre-period distribution, and
returns the post-period value at that same rank within group $g$’s
distribution. The assumption requires this operation to give the same
result for both groups at every $y$. It can be thought of as a nonlinear
version of parallel trends: where parallel trends requires the *mean* of
$Y(0)$ to shift equally across groups, CiC requires the *entire
distributional time-change* to be the same.

**Identification:** Under this assumption the counterfactual quantile
function is
$$Q_{Y_{t}(0)|D = 1}(\tau) = Q_{Y_{t}|D = 0}\!\left( F_{Y_{t - 1}|D = 0}\!\left( Q_{Y_{t - 1}|D = 1}(\tau) \right) \right).$$

**When to use:** Works with panel or repeated cross sections.

``` r
res_cic <- cic(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 100
)
summary(res_cic)
```

    #>
    #> Overall ATT:
    #>      ATT    Std. Error     [ 95%  Conf. Int.]
    #>  -0.0197        0.0157    -0.0541      0.0147
    #>
    #>
    #> Dynamic Effects:
    #>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]
    #>          -3   0.0508     0.0214        0.0090      0.0927 *
    #>          -2   0.0158     0.0169       -0.0174      0.0491
    #>          -1  -0.0128     0.0179       -0.0478      0.0222
    #>           0  -0.0081     0.0154       -0.0383      0.0221
    #>           1  -0.0364     0.0233       -0.0820      0.0092
    #>           2  -0.1226     0.0500       -0.2206     -0.0245 *
    #>           3  -0.0930     0.0502       -0.1914      0.0055
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

------------------------------------------------------------------------

## Quantile Difference-in-Differences — `qdid()`

**References:** Athey and Imbens ([2006](#ref-athey-imbens-2006)); Meyer
et al. ([1995](#ref-meyer-viscusi-durbin-1995))

**Identifying assumption:** The change in quantiles over time is the
same for both groups:
$$Q_{Y_{t}(0)|D = 1}(\tau) - Q_{Y_{t - 1}(0)|D = 1}(\tau) = Q_{Y_{t}(0)|D = 0}(\tau) - Q_{Y_{t - 1}(0)|D = 0}(\tau)\quad{\text{for all}\mspace{6mu}}\tau \in (0,1).$$

**Intuition:** This requires parallel trends to hold at every quantile
level $\tau$, not just in the mean. It directly identifies the
counterfactual:
$$Q_{Y_{t}(0)|D = 1}(\tau) = Q_{Y_{t - 1}|D = 1}(\tau) + Q_{Y_{t}|D = 0}(\tau) - Q_{Y_{t - 1}|D = 0}(\tau).$$

**When to use:** Works with panel or repeated cross sections.

``` r
res_qdid <- qdid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 100
)
summary(res_qdid)
```

    #>
    #> Overall ATT:
    #>      ATT    Std. Error     [ 95%  Conf. Int.]
    #>  -0.0282        0.0125    -0.0506     -0.0058 *
    #>
    #>
    #> Dynamic Effects:
    #>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]
    #>          -3   0.0266     0.0172       -0.0070      0.0602
    #>          -2  -0.0040     0.0133       -0.0301      0.0221
    #>          -1  -0.0247     0.0126       -0.0493     -0.0001 *
    #>           0  -0.0164     0.0128       -0.0415      0.0086
    #>           1  -0.0584     0.0189       -0.0953     -0.0214 *
    #>           2  -0.1470     0.0408       -0.2269     -0.0671 *
    #>           3  -0.1111     0.0404       -0.1904     -0.0319 *
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

------------------------------------------------------------------------

## Panel QTT — `panel_qtt()`

**Reference:** Callaway and Li ([2019](#ref-callaway-li-2019))

**Identifying assumptions:** The first assumption is
$\Delta Y_{t}(0) \perp D$ — the distribution of the change in untreated
potential outcomes is the same for both groups:
$$F_{\Delta Y_{t}(0)|D = 1} = F_{\Delta Y_{t}(0)|D = 0}.$$ This
identifies the marginal distribution of $\Delta Y_{t}(0)$ for the
treated group. As shown in Fan and Yu ([2012](#ref-fan-yu-2012)), this
is not sufficient on its own to identify the QTT. Identification is
completed by a copula stability assumption: the copula linking
$\Delta Y_{t}(0)$ to the lagged outcome $Y_{t - 1}(0)$ is stable over
time within the treated group,
$$C_{\Delta Y_{t}(0),\, Y_{t - 1}(0)|D = 1} = C_{\Delta Y_{t - 1}(0),\, Y_{t - 2}(0)|D = 1}.$$

**Intuition:** For treated units, $Y_{t - 1}$ and $Y_{t - 2}$ are
observed pre-treatment, so the right-hand copula can be estimated
directly from the data. Copula stability says this same dependence
structure governs $(\Delta Y_{t}(0),Y_{t - 1}(0))$ in the post-period.
Combined with the identified marginals, this pins down the conditional
distribution of $\Delta Y_{t}(0)$ given $Y_{t - 1}(0) = y$ for treated
units. Since $Y_{t}(0) = Y_{t - 1} + \Delta Y_{t}(0)$ and $Y_{t - 1}$ is
observed, the QTT is identified. Requires panel data with at least two
pre-treatment periods per cohort.

``` r
res_pqtt <- panel_qtt(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 100
)
summary(res_pqtt)
```

    #>
    #> Overall QTT Curve:
    #>  Quantile     QTT Std. Error 95% CB Lower 95% CB Upper
    #>       0.1 -0.0270     0.1625      -0.6175       0.5634
    #>       0.2 -0.0852     0.0833      -0.3879       0.2176
    #>       0.3 -0.0250     0.0523      -0.2149       0.1650
    #>       0.4 -0.0680     0.0612      -0.2905       0.1546
    #>       0.5 -0.0965     0.0545      -0.2944       0.1014
    #>       0.6  0.0070     0.0375      -0.1294       0.1433
    #>       0.7  0.0409     0.0567      -0.1652       0.2470
    #>       0.8 -0.0201     0.0581      -0.2312       0.1910
    #>       0.9 -0.0832     0.0632      -0.3130       0.1466

------------------------------------------------------------------------

## Distributional DiD — `ddid()`

**Reference:** Callaway and Li ([2019](#ref-callaway-li-2019))

**Identifying assumptions:** The first is $\Delta Y_{t}(0) \perp D$ —
the same as in
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md).
The copula assumption differs: rather than stability over time,
[`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) assumes
the copula of $(\Delta Y_{t}(0),Y_{t - 1}(0))$ is the same for treated
and untreated groups,
$$C_{\Delta Y_{t}(0),\, Y_{t - 1}(0)|D = 1} = C_{\Delta Y_{t}(0),\, Y_{t - 1}(0)|D = 0}.$$

**Intuition:** For the untreated group, $\Delta Y_{t}(0) = \Delta Y_{t}$
is observed, so the right-hand copula can be estimated directly. The
assumption transfers it to the treated group. Combined with the
identified marginals, this pins down the joint distribution of
$(\Delta Y_{t}(0),Y_{t - 1}(0))$ for treated units and identifies the
QTT — in exactly the same way as
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md).
The key contrast is the source of the copula:
[`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) borrows
it from the untreated group at the same point in time, while
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
borrows it from the treated group’s own pre-treatment history.

``` r
res_ddid <- ddid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 100
)
summary(res_ddid)
```

    #>
    #> Overall ATT:
    #>      ATT    Std. Error     [ 95%  Conf. Int.]
    #>  -0.0453        0.0143    -0.0785      -0.012 *
    #>
    #>
    #> Dynamic Effects:
    #>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]
    #>          -3   0.0134     0.0159       -0.0178      0.0446
    #>          -2  -0.0146     0.0151       -0.0442      0.0150
    #>          -1  -0.0403     0.0154       -0.0704     -0.0102 *
    #>           0  -0.0324     0.0132       -0.0583     -0.0065 *
    #>           1  -0.0637     0.0182       -0.0994     -0.0281 *
    #>           2  -0.1377     0.0352       -0.2068     -0.0686 *
    #>           3  -0.1087     0.0367       -0.1806     -0.0368 *
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

------------------------------------------------------------------------

## Mean Difference-in-Differences — `mdid()`

**Reference:** Thuysbaert ([2007](#ref-thuysbaert-2007))

**Identifying assumption:** The counterfactual quantile for the treated
group in the post-period is the pre-period quantile shifted by the mean
change in outcomes for the untreated group:
$$Q_{Y_{t}(0)|D = 1}(\tau) = Q_{Y_{t - 1}|D = 1}(\tau) + E\lbrack Y_{t}|D = 0\rbrack - E\lbrack Y_{t - 1}|D = 0\rbrack.$$

**Intuition:** The assumption says the shape of the treated group’s
outcome distribution is unchanged in the counterfactual — the entire
distribution simply shifts by the amount the untreated group’s mean
changed. This implies mean parallel trends (so the ATT is identified),
but is stronger: it additionally requires that only the location of the
distribution shifts, not its shape. Works with panel or repeated cross
sections.

``` r
res_mdid <- mdid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 100
)
summary(res_mdid)
```

    #>
    #> Overall ATT:
    #>      ATT    Std. Error     [ 95%  Conf. Int.]
    #>  -0.0305        0.0114    -0.0496     -0.0113 *
    #>
    #>
    #> Dynamic Effects:
    #>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]
    #>          -3   0.0298     0.0127        0.0048      0.0547 *
    #>          -2  -0.0024     0.0158       -0.0334      0.0285
    #>          -1  -0.0243     0.0131       -0.0500      0.0014
    #>           0  -0.0189     0.0106       -0.0397      0.0018
    #>           1  -0.0536     0.0173       -0.0875     -0.0197 *
    #>           2  -0.1363     0.0393       -0.2133     -0.0593 *
    #>           3  -0.1008     0.0378       -0.1748     -0.0268 *
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

------------------------------------------------------------------------

## Lagged-Outcome Unconfoundedness — `lou_qtt()`

**Identification assumption:** Unconfoundedness conditional on the
lagged outcome — treatment is as good as randomly assigned after
conditioning on the pre-treatment value of the outcome.

**When to use:** When selection into treatment may depend on lagged
outcomes but not on unobserved heterogeneity — i.e., once you condition
on the pre-treatment outcome, treatment assignment is as good as random.
Requires panel data.

``` r
res_lou <- lou_qtt(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "att", biters = 100
)
summary(res_lou)
```

    #>
    #> Overall ATT:
    #>      ATT    Std. Error     [ 95%  Conf. Int.]
    #>  -0.0384        0.0166    -0.0725     -0.0043 *
    #>
    #>
    #> Dynamic Effects:
    #>  Event Time Estimate Std. Error [95% Simult.  Conf. Band]
    #>          -3   0.0294     0.0144        0.0012      0.0577 *
    #>          -2  -0.0124     0.0161       -0.0441      0.0192
    #>          -1  -0.0310     0.0199       -0.0701      0.0081
    #>           0  -0.0243     0.0157       -0.0551      0.0065
    #>           1  -0.0869     0.0306       -0.1470     -0.0269 *
    #>           2  -0.1394     0.0356       -0.2092     -0.0696 *
    #>           3  -0.1093     0.0378       -0.1834     -0.0352 *
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

------------------------------------------------------------------------

## Summary

| Function                                                                  | Assumption                      | Panel required | Notes                               |
|---------------------------------------------------------------------------|---------------------------------|----------------|-------------------------------------|
| [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md)             | Rank invariance + support       | Optional       | Nonparametric; strongest assumption |
| [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md)           | Distributional parallel trends  | Optional       | DiD analog for distributions        |
| [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md) | Copula stability                | Yes            | Requires ≥ 2 pre-periods per cohort |
| [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md)           | Conditional distributional PT   | Yes            | Adds covariate adjustment to QDiD   |
| [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md)           | Mean parallel trends            | Optional       | Weakest DiD assumption              |
| [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)     | Lagged-outcome unconfoundedness | Yes            | Not a parallel trends assumption    |

All estimators return either a `pte_results` object (`gt_type = "att"`)
or a `pte_qtt` object (`gt_type = "qtt"`), and support
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
for plotting overall curves and event-study summaries. See
[`vignette("staggered-adoption")`](https://bcallaway11.github.io/qte/articles/staggered-adoption.md)
for the full applied workflow.

## References

Athey, Susan, and Guido Imbens. 2006. “Identification and Inference in
Nonlinear Difference-in-Differences Models.” *Econometrica* 74 (2):
431–97.

Callaway, Brantly, and Tong Li. 2019. “Quantile Treatment Effects in
Difference in Differences Models with Panel Data.” *Quantitative
Economics* 10 (4): 1579–618.

Fan, Yanqin, and Zhengfei Yu. 2012. “Partial Identification of
Distributional and Quantile Treatment Effects in
Difference-in-Differences Models.” *Economics Letters* 115 (3): 511–15.

Meyer, Bruce, Kip Viscusi, and David Durbin. 1995. “Workers’
Compensation and Injury Duration: Evidence from a Natural Experiment.”
*American Economic Review*, 322–40.

Thuysbaert, Bram. 2007. “Distributional Comparisons in Difference in
Differences Models.” Unpublished manuscript.
