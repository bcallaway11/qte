# Panel Data Estimators for Quantile Treatment Effects

## Overview

The `qte` package provides seven estimators of the quantile treatment
effect on the treated (QTT) for settings where units are observed before
and after treatment. All share the same `yname`/`gname`/`tname`/`idname`
interface and support staggered treatment adoption. They differ in their
identifying assumptions. Six are difference-in-differences-type
assumptions, which restrict how untreated outcomes change over time;
one,
[`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md),
is a selection-on-observables assumption that conditions on the lagged
outcome.

This vignette states each assumption and shows the minimal call. For the
full applied workflow (aggregations, event studies, testing whether the
effect is more than a location shift, covariates), see
[`vignette("staggered-adoption")`](https://bcallaway11.github.io/qte/articles/staggered-adoption.md).

``` r
library(qte)
library(ptetools)
library(ggplot2)
set.seed(42)
data(mpdta, package = "did")
```

We use `mpdta` throughout: a balanced panel of 500 US counties observed
annually from 2003–2007, with `lemp` (log employment) as the outcome and
`first.treat` as the cohort variable (year of first minimum wage
increase; 0 = never treated).

The common call structure for all seven estimators is:

``` r
estimator(
  yname  = "lemp",        # outcome
  gname  = "first.treat", # first treatment period (0 = never treated)
  tname  = "year",        # calendar time
  idname = "countyreal",  # unit id
  data   = mpdta,
  gt_type = "qtt",        # "qtt" (QTT curve) or "att"
  probs   = seq(0.1, 0.9, 0.1),
  biters  = 100
)
```

With `gt_type = "qtt"`, each estimator returns the QTT curve together
with the ATT implied by the same assumption, both from the same
bootstrap draws. Each plot below shows the overall QTT curve with a
uniform confidence band and the overall ATT as a horizontal line.

All of these estimators except
[`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
are built on quantiles or ranks and are designed for continuous
outcomes. With a discrete outcome, or one with a large mass point (such
as many zeros), quantile treatment effects are usually not point
identified, and the estimators warn.

------------------------------------------------------------------------

## Change in Changes — `cic()`

**Reference:** Athey and Imbens ([2006](#ref-athey-imbens-2006))

**Identifying assumption:** The distributional time-change is the same
for both groups:
``` math
Q_{Y_t(0)|D=1}\!\left(F_{Y_{t-1}(0)|D=1}(y)\right) = Q_{Y_t(0)|D=0}\!\left(F_{Y_{t-1}(0)|D=0}(y)\right) \quad \text{for all } y,
```
along with a support condition:
$`\text{supp}(Y_{t-1}(0)|D=1) \subseteq \text{supp}(Y_{t-1}(0)|D=0)`$.

**Intuition:** The expression $`Q_{Y_t(0)|D=g}(F_{Y_{t-1}(0)|D=g}(y))`$
takes a pre-period value $`y`$, finds its quantile rank in group $`g`$’s
pre-period distribution, and returns the post-period value at that same
rank within group $`g`$’s distribution. The assumption requires this
operation to give the same result for both groups at every $`y`$. It can
be thought of as a nonlinear version of parallel trends: where parallel
trends requires the *mean* of $`Y(0)`$ to shift equally across groups,
CiC requires the *entire distributional time-change* to be the same.
Athey and Imbens ([2006](#ref-athey-imbens-2006)) derive it from a model
$`Y(0) = h(U, t)`$ with $`h`$ increasing in unobservables $`U`$ whose
distribution does not change over time within a group.

**Identification:** Under this assumption the counterfactual quantile
function is
``` math
Q_{Y_t(0)|D=1}(\tau) = Q_{Y_t|D=0}\!\left(F_{Y_{t-1}|D=0}\!\left(Q_{Y_{t-1}|D=1}(\tau)\right)\right).
```

**Data:** Panel data or repeated cross sections (`panel = FALSE`).

``` r
res_cic <- cic(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 100
)
autoplot(res_cic)
```

![Overall QTT curve from
cic()](panel-estimators_files/figure-html/cic-plot-1.png)

------------------------------------------------------------------------

## Quantile Difference-in-Differences — `qdid()`

**References:** Athey and Imbens ([2006](#ref-athey-imbens-2006));
Meyer, Viscusi, and Durbin ([1995](#ref-meyer-viscusi-durbin-1995))

**Identifying assumption:** The change in quantiles over time is the
same for both groups:
``` math
Q_{Y_t(0)|D=1}(\tau) - Q_{Y_{t-1}(0)|D=1}(\tau) = Q_{Y_t(0)|D=0}(\tau) - Q_{Y_{t-1}(0)|D=0}(\tau) \quad \text{for all } \tau \in (0,1).
```

**Intuition:** This requires parallel trends to hold at every quantile
level $`\tau`$, not just in the mean. It directly identifies the
counterfactual:
``` math
Q_{Y_t(0)|D=1}(\tau) = Q_{Y_{t-1}|D=1}(\tau) + Q_{Y_t|D=0}(\tau) - Q_{Y_{t-1}|D=0}(\tau).
```

**Data:** Panel data or repeated cross sections (`panel = FALSE`).

``` r
res_qdid <- qdid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 100
)
autoplot(res_qdid)
```

![Overall QTT curve from
qdid()](panel-estimators_files/figure-html/qdid-plot-1.png)

------------------------------------------------------------------------

## Panel QTT — `panel_qtt()`

**Reference:** Callaway and Li ([2019](#ref-callaway-li-2019))

**Identifying assumptions:** The first assumption is
$`\Delta Y_t(0) \perp D`$ — the distribution of the change in untreated
potential outcomes is the same for both groups:
``` math
F_{\Delta Y_t(0)|D=1} = F_{\Delta Y_t(0)|D=0}.
```
This identifies the marginal distribution of $`\Delta Y_t(0)`$ for the
treated group. As shown in Fan and Yu ([2012](#ref-fan-yu-2012)), this
is not sufficient on its own to identify the QTT. Identification is
completed by a copula stability assumption: the copula linking
$`\Delta Y_t(0)`$ to the lagged outcome $`Y_{t-1}(0)`$ is stable over
time within the treated group,
``` math
C_{\Delta Y_t(0),\, Y_{t-1}(0)|D=1} = C_{\Delta Y_{t-1}(0),\, Y_{t-2}(0)|D=1}.
```

**Intuition:** For treated units, $`Y_{t-1}`$ and $`Y_{t-2}`$ are
observed pre-treatment, so the right-hand copula can be estimated
directly from the data. Copula stability says this same dependence
structure governs $`(\Delta Y_t(0), Y_{t-1}(0))`$ in the post-period.
Combined with the identified marginals, this pins down the conditional
distribution of $`\Delta Y_t(0)`$ given $`Y_{t-1}(0) = y`$ for treated
units. Since $`Y_t(0) = Y_{t-1} + \Delta Y_t(0)`$ and $`Y_{t-1}`$ is
observed, the QTT is identified.

**Data:** Panel data, with at least two pre-treatment periods per
cohort. With staggered adoption, the default `pre_copula = "long"` uses
a pre-treatment window as long as the post-treatment horizon, so longer
event times need more pre-treatment periods; `pre_copula = "short"`
always uses the two periods just before treatment.

``` r
res_pqtt <- panel_qtt(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 100
)
autoplot(res_pqtt)
```

![Overall QTT curve from
panel_qtt()](panel-estimators_files/figure-html/panel-qtt-plot-1.png)

------------------------------------------------------------------------

## Distributional DiD — `ddid()`

**Reference:** Callaway and Li ([2019](#ref-callaway-li-2019))

**Identifying assumptions:** The first is $`\Delta Y_t(0) \perp D`$ —
the same as in
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md).
The copula assumption differs: rather than stability over time,
[`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) assumes
the copula of $`(\Delta Y_t(0), Y_{t-1}(0))`$ is the same for treated
and untreated groups,
``` math
C_{\Delta Y_t(0),\, Y_{t-1}(0)|D=1} = C_{\Delta Y_t(0),\, Y_{t-1}(0)|D=0}.
```

**Intuition:** For the untreated group, $`\Delta Y_t(0) = \Delta Y_t`$
is observed, so the right-hand copula can be estimated directly. The
assumption transfers it to the treated group. Combined with the
identified marginals, this pins down the joint distribution of
$`(\Delta Y_t(0), Y_{t-1}(0))`$ for treated units and identifies the QTT
— in exactly the same way as
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md).
The key contrast is the source of the copula:
[`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) borrows
it from the untreated group at the same point in time, while
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
borrows it from the treated group’s own pre-treatment history.

**Data:** Panel data, with one pre-treatment period per cohort.

``` r
res_ddid <- ddid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 100
)
autoplot(res_ddid)
```

![Overall QTT curve from
ddid()](panel-estimators_files/figure-html/ddid-plot-1.png)

------------------------------------------------------------------------

## Mean Difference-in-Differences — `mdid()`

**Reference:** Thuysbaert ([2007](#ref-thuysbaert-2007))

**Identifying assumption:** The counterfactual quantile for the treated
group in the post-period is the pre-period quantile shifted by the mean
change in outcomes for the untreated group:
``` math
Q_{Y_t(0)|D=1}(\tau) = Q_{Y_{t-1}|D=1}(\tau) + E[Y_t|D=0] - E[Y_{t-1}|D=0].
```

**Intuition:** The assumption says the shape of the treated group’s
outcome distribution is unchanged in the counterfactual — the entire
distribution simply shifts by the amount the untreated group’s mean
changed. This implies mean parallel trends, so the ATT is the usual DiD
estimate, but it is stronger: it additionally requires that only the
location of the untreated distribution changes over time, not its shape.

**Data:** Panel data or repeated cross sections (`panel = FALSE`).

``` r
res_mdid <- mdid(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 100
)
autoplot(res_mdid)
```

![Overall QTT curve from
mdid()](panel-estimators_files/figure-html/mdid-plot-1.png)

------------------------------------------------------------------------

## Parallel Trends in Distribution Functions — `cdf_did()`

**Reference:** Roth and Sant’Anna ([2023](#ref-roth-santanna-2021))

**Identifying assumption:** The change over time in the distribution
function of untreated potential outcomes is the same for both groups, at
every value $`y`$:
``` math
F_{Y_t(0)|D=1}(y) - F_{Y_{t-1}(0)|D=1}(y) = F_{Y_t(0)|D=0}(y) - F_{Y_{t-1}(0)|D=0}(y) \quad \text{for all } y.
```

**Intuition:** This is parallel trends for every indicator
$`1\{Y \le y\}`$ instead of for $`Y`$ itself. Roth and Sant’Anna
([2023](#ref-roth-santanna-2021)) show that it is essentially the
condition under which parallel trends holds for every strictly monotone
transformation of the outcome (logs, levels, …), so it is the natural
assumption when the choice of units is arbitrary. It implies mean
parallel trends, so the ATT is the usual DiD estimate. The
counterfactual distribution is
``` math
F_{Y_t(0)|D=1}(y) = F_{Y_{t-1}|D=1}(y) + F_{Y_t|D=0}(y) - F_{Y_{t-1}|D=0}(y).
```
In finite samples the right-hand side need not be a valid distribution
function;
[`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
rearranges it into one ([Chernozhukov, Fernandez-Val, and Galichon
2010](#ref-chernozhukov-val-galichon-2010)) and reports how much that
changed it.

Unlike the rank- and quantile-based estimators above, this assumption
identifies the counterfactual distribution for discrete outcomes too.
With covariates, `est_method = c("aipw", "ipw", "reg")` gives doubly
robust, inverse-probability-weighted and outcome-regression versions of
the conditional assumption ([Sant’Anna and Zhao
2020](#ref-santanna-zhao-2020)), applied at each threshold $`y`$.

**Data:** Panel data or repeated cross sections (`panel = FALSE`).

``` r
res_cdf <- cdf_did(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 100
)
autoplot(res_cdf)
```

![Overall QTT curve from
cdf_did()](panel-estimators_files/figure-html/cdf-did-plot-1.png)

------------------------------------------------------------------------

## Lagged-Outcome Unconfoundedness — `lou_qtt()`

**Identifying assumption:** Conditional on the pre-treatment outcome
(and any covariates $`X`$), untreated potential outcomes are independent
of treatment:
``` math
Y_t(0) \perp D \mid Y_{t-1}, X.
```

**Intuition:** Treatment is as good as randomly assigned among units
with the same pre-treatment outcome. This is not a parallel trends
assumption: selection may depend on the level of the lagged outcome but,
given it, not on unobserved unit-specific heterogeneity. DiD-type
assumptions instead allow selection on time-invariant unobservables. The
two are not nested, and when they disagree it is informative to report
both. The counterfactual distribution is estimated by inverse
probability weighting, outcome (quantile) regression, or a doubly robust
combination (`est_method = c("aipw", "ipw", "reg")`).

**Data:** Panel data. With `lagged_outcome_cov = FALSE` the lagged
outcome is dropped from the conditioning set, leaving unconfoundedness
given $`X`$, which also works with repeated cross sections.

``` r
res_lou <- lou_qtt(
  yname  = "lemp", gname = "first.treat", tname = "year",
  idname = "countyreal", data = mpdta,
  gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 100
)
autoplot(res_lou)
```

![Overall QTT curve from
lou_qtt()](panel-estimators_files/figure-html/lou-plot-1.png)

------------------------------------------------------------------------

## Summary

| Function | Identifying assumption | Repeated cross sections | Notes |
|----|----|----|----|
| [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) | Same distributional time-change for both groups (changes-in-changes) | Yes | Needs a support condition; continuous outcomes |
| [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md) | Parallel trends in every quantile | Yes | Continuous outcomes |
| [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md) | Parallel trends in the distribution of changes + copula stable over time | No | Needs ≥ 2 pre-treatment periods |
| [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) | Parallel trends in the distribution of changes + same copula for both groups | No |  |
| [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) | Treated distribution shifts by the comparison group’s mean change | Yes | Implies mean parallel trends, and is stronger |
| [`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md) | Parallel trends in distribution functions | Yes | Also identified for discrete outcomes |
| [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md) | Unconfoundedness given the lagged outcome | With `lagged_outcome_cov = FALSE` | Not a parallel trends assumption |

All seven take covariates through `xformula`. All return a `pte_qtt`
object for `gt_type = "qtt"`, or a `pte_emp_boot` object for
`gt_type = "att"`, with
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods. See
[`vignette("staggered-adoption")`](https://bcallaway11.github.io/qte/articles/staggered-adoption.md)
for the full applied workflow.

## References

Athey, Susan, and Guido Imbens. 2006. “Identification and Inference in
Nonlinear Difference-in-Differences Models.” *Econometrica* 74 (2):
431–97.

Callaway, Brantly, and Tong Li. 2019. “Quantile Treatment Effects in
Difference in Differences Models with Panel Data.” *Quantitative
Economics* 10 (4): 1579–1618.

Chernozhukov, Victor, Ivan Fernandez-Val, and Alfred Galichon. 2010.
“Quantile and Probability Curves Without Crossing.” *Econometrica* 78
(3): 1093–1125.

Fan, Yanqin, and Zhengfei Yu. 2012. “Partial Identification of
Distributional and Quantile Treatment Effects in
Difference-in-Differences Models.” *Economics Letters* 115 (3): 511–15.

Meyer, Bruce, Kip Viscusi, and David Durbin. 1995. “Workers’
Compensation and Injury Duration: Evidence from a Natural Experiment.”
*American Economic Review*, 322–40.

Roth, Jonathan, and Pedro HC Sant’Anna. 2023. “When Is Parallel Trends
Sensitive to Functional Form?” *Econometrica* 91 (2): 737–47.

Sant’Anna, Pedro HC, and Jun Zhao. 2020. “Doubly Robust
Difference-in-Differences Estimators.” *Journal of Econometrics* 219
(1): 101–22.

Thuysbaert, Bram. 2007. “Distributional Comparisons in Difference in
Differences Models.”
