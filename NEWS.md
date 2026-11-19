# qte 2.1.0 (development version)

## Breaking changes

* Removed the deprecated two-period functions, as announced in 2.0.0:
  `CiC()`, `QDiD()`, `MDiD()`, `ddid2()`, `panel.qtet()`, `ci.qte()`,
  `ci.qtet()` and `ggqte()`. Use `cic()`, `qdid()`, `mdid()`, `ddid()`,
  `panel_qtt()`, `unc_qte()` (with `target = "qtt"` for the QTT; `lou_qtt()`
  for staggered adoption) and `autoplot()` instead.
* Removed the exported helper `panelize.data()`, which only served those
  functions.
* `formula.tools` is no longer a dependency.

## New features

* The covariate-adjusted estimators gain `cond_dist = c("qr", "dr")`:
  `cic()`, `qdid()`, `ddid()`, `panel_qtt()` (and their `*_gt()` workhorses),
  and the outcome-regression and AIPW branches of `lou_qtt()` and `unc_qte()`.
  With covariates, `"dr"` estimates conditional outcome distributions by
  distribution regression (a logit of 1{Y <= y} on the covariates at a grid of
  thresholds) instead of quantile regression. It handles outcomes with mass
  points (e.g. many zeros) naturally and avoids truncating each unit's
  conditional distribution at its 1st and 99th percentiles. The default stays
  `"qr"`, and the `"qr"` results are unchanged.

## Bug fixes (these change results)

* `ddid()` / `ddid_gt()` with covariates now integrate the counterfactual
  over the *treated* group's covariate distribution, by reweighting control
  units with the propensity-score odds p(X)/(1 - p(X)) (#10). Previously the
  controls' covariate distribution was used, which biased the ATT and QTT
  whenever treatment depends on the covariates (bias 0.84 in a simple
  simulated design; now about 0). Results without covariates are unchanged.
* `est_method = "aipw"` in `unc_qte()` and `lou_qtt()` now reports a
  doubly-robust ATT/ATE: the mean of the AIPW counterfactual distribution,
  using the outcome model's conditional mean (#11). Previously the reported
  ATT/ATE was the IPW estimate, so it ignored the outcome model and wasn't
  doubly robust. The AIPW quantile curves are unchanged.
* `unc_qte()` now uses `weightsname` throughout: in the propensity-score
  model, in the outcome-regression fits, and in the AIPW distribution
  functions (previously only in the final averaging). `lou_qtt()`'s
  propensity-score model is now weighted too. With integer weights, the
  estimates now equal those from the unweighted fit on row-replicated data.
  Results without weights are unchanged.

* With groups and time periods coded as calendar years, naming the group
  column `G` or the time column `period` gave `NaN` estimates (ptetools
  #13). Fixed in ptetools, which now keeps its internal columns separate
  from the user's; `panel_qtt()`'s subset function uses those internal
  columns. Requires `ptetools (>= 1.0.2)`. Results for other column names
  are unchanged.

* `panel_qtt()` now treats `anticipation` in units of the time variable, as
  `did::att_gt()` and ptetools do: a cell is post-treatment from
  `g - anticipation` on, pre1 is the last period at least `anticipation`
  before `g`, and not-yet-treated comparison units must not be treated or
  anticipating treatment by the cell's latest period. It previously counted
  anticipation in periods and ignored it for the comparison group. Results
  without anticipation are unchanged.

* `panel_qtt(gt_type = "qtt")` with the default `pre_copula = "long"` no
  longer fails ("no valid post-treatment cells") or drops valid cells when
  time periods are unequally spaced, e.g. `lalonde.psid.panel` (1974, 1975,
  1978). The aggregation step checked whether the earlier base period exists
  in calendar years, while the estimator chooses base periods by position in
  the panel. Results with equally spaced periods are unchanged.

## Minor improvements

* `unc_qte()` now checks up front that `dname` is a column of `data` coded
  0/1, with an informative error. A treatment coded e.g. 1/2 previously
  failed inside the propensity-score fit with
  "y values must be 0 <= y <= 1" (#1).

# qte 2.0.1

* Fixed covariate-adjusted `cic()`, `qdid()`, `ddid()` and `panel_qtt()`
  (and their `*_gt()` workhorses) (#9):
  * They no longer fail when a group has a single treated unit.
  * Conditional quantiles are now rearranged (sorted) before use, so the
    conditional CDFs are monotone when quantile regressions cross.
  * Conditional ranks are no longer one quantile-grid step (0.01) too low.
    This came from quantreg's `predict(..., type = "Fhat", stepfun = TRUE)`
    and biased the counterfactual distribution downward.
  * Conditional CDFs and quantile functions are now linearly interpolated
    between the estimated quantiles, not step functions, so ranks are no
    longer rounded down to the quantile grid. Beyond the first and last
    estimated quantiles they stay flat (no tail extrapolation).

  Covariate-adjusted estimates change as a result. In a simulation where
  conditional CiC holds exactly, the bias of the ATT falls from 0.15 to 0.02.
  Estimates without covariates are unchanged.

* `cic()`, `qdid()`, `ddid()`, `mdid()`, `panel_qtt()` and `lou_qtt()` are
  faster, with the same results:
  * Weighted ranks and quantiles in the `*_gt()` functions use one sort
    instead of a pass over the data for each unit, so they are no longer
    quadratic in the group size.
  * The covariate versions of `cic()` and `qdid()` fit the post-period
    quantile regression for the comparison group once per period, not once
    per cohort, because it is the same for every cohort already treated by
    that period.
  * Only the columns the estimators use are passed to `ptetools::pte()`.

  Combined with the faster empirical bootstrap in ptetools 1.0.2,
  unconditional CiC on a 3,000-unit, 30-period panel runs more than 20 times
  faster.

# qte 2.0.0

This is a major update that extends all estimators to support staggered
treatment adoption. The API has changed substantially — please read this
section before upgrading.

## New functions

* `unc_qte()` — cross-sectional QTE and QTT under unconfoundedness. Supports
  IPW, outcome regression, and doubly robust (AIPW) estimation. Uses a
  `yname`/`dname` interface consistent with the `did`/`ptetools` convention.
  Replaces `ci.qte()` and `ci.qtet()`, which are retained as deprecated
  wrappers.

* `cic()` / `cic_gt()` — Change in Changes estimator for staggered treatment
  adoption (Athey and Imbens 2006). Supports panel and repeated cross sections
  data, and returns either ATT (`gt_type = "att"`) or a full QTT curve
  (`gt_type = "qtt"`). Replaces `CiC()`, which is retained as a deprecated
  wrapper.

* `qdid()` / `qdid_gt()` — Quantile Difference in Differences estimator for
  staggered treatment adoption (Callaway, Li, and Oka 2018). Replaces `QDiD()`,
  which is retained as a deprecated wrapper.

* `mdid()` / `mdid_gt()` — Mean Difference in Differences estimator for
  staggered treatment adoption (Callaway and Li 2019). Replaces `MDiD()`,
  which is retained as a deprecated wrapper.

* `ddid()` / `ddid_gt()` — Distributional Difference in Differences estimator
  for staggered treatment adoption (Callaway and Li 2019). Replaces `ddid2()`,
  which is retained as a deprecated wrapper.

* `panel_qtt()` / `panel_qtt_gt()` — Panel QTT estimator using copula transfer
  across pre-treatment periods (Callaway and Li 2019), with staggered treatment
  adoption support. The `pre_copula` argument controls whether a rolling
  (`"long"`, default) or fixed (`"short"`) base period is used to transfer the
  copula. Replaces `panel.qtet()`, which is retained as a deprecated wrapper.

* `lou_qtt()` — Lagged-outcome unconfoundedness QTT estimator for staggered
  treatment adoption. Uses the lagged outcome as a control variable to identify
  the QTT. New; no prior equivalent in the package.

## Deprecated functions

The following functions emit a runtime deprecation warning via `.Deprecated()`
and will be removed in v2.1 (target: 2026-08). Each points to its modern
replacement in the warning message.

| Deprecated | Replacement |
|---|---|
| `ci.qte()` | `unc_qte(target = "qte")` |
| `ci.qtet()` | `unc_qte(target = "qtt")` |
| `CiC()` | `cic()` |
| `QDiD()` | `qdid()` |
| `MDiD()` | `mdid()` |
| `ddid2()` | `ddid()` |
| `panel.qtet()` | `panel_qtt()` |

## Removed functions

* `computeDiffSE()` — was dead code; never called anywhere in the package.
  Design notes for a future QTT-vs-ATT test are in `dev/NOTES.md`.

* `DR()`, `spatt()`, `bounds()` — removed in earlier 1.5.x development work.

* `qtes2mat()`, `diffQ()` — formatting helpers tied to the legacy `QTE` class
  with no internal callers.

## Infrastructure

* Full `testthat` test suite added (83 fast tests; 10 slow bootstrap tests
  skipped by default, enabled via `R_SLOW_TESTS=true`). The package previously
  had no tests.

* All new estimators are backed by `ptetools::pte()`. The `attgt_fun` hook
  pattern makes it straightforward to add new identification strategies without
  duplicating bootstrap or aggregation code.

* Estimators now pass an explicit `aggte_fun` to `ptetools::pte()`,
  replacing the earlier implicit dispatch based on `gt_type`.

* `ptetools` added to `Imports` (was missing from `DESCRIPTION` despite being
  used).

* `Hmisc` and `texreg` removed from `Imports` (were unused after earlier
  cleanup).

* License updated from GPL-2 to GPL-3.

---

# qte 1.4.0

* Added `cic2()` for Change in Changes with multiple periods and variation in
  treatment timing, using `ptetools` as a backend. This function was
  experimental; it has been superseded by `cic()` in v2.0.0.

# qte 1.3.3

* Bug fix: critical values for uniform confidence bands now correctly reflect
  the `alp` argument rather than being stuck at `alp = 0.05`.

# qte 1.3.2

* Significance level is now passed through in the results.
* Added ability to change significance level in `ggqte()` for pointwise
  confidence intervals without re-estimating.

# qte 1.3.1

* Updated maintainer contact information.
* Improvements and updates to vignettes.

# qte 1.3.0

* Added first-step quantile regression support for `panel.qtet()`, `QDiD()`,
  and `CiC()`.
* Added checks that input data is in the correct format, particularly for panel
  data.

# qte 1.2.3

* Now using empirical quantiles in the CiC method.
* Added code for computing uniform confidence bands.
* General improvements to the codebase, including combining shared code across
  methods.

# qte 1.2.2

* Bug fixes to the CiC method.

# qte 1.2

* Added `ddid2()` to compute quantile treatment effects under a Difference in
  Differences assumption following Callaway and Li (2016).
* Covariates for all methods now specified using formulas (`xformla` argument)
  rather than passing variable names directly.
* Improvements to bootstrapping standard errors.
* Some extensions to computing the ATT via `spatt()`.
* Preliminary version of multi-period ATT moved to the `did` package.
