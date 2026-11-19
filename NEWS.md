# qte 2.1.0 (development version)

## Breaking changes

* Removed the deprecated two-period functions, as announced in 2.0.0:
  `CiC()`, `QDiD()`, `MDiD()`, `ddid2()`, `panel.qtet()`, `ci.qte()`,
  `ci.qtet()` and `ggqte()`. Use `cic()`, `qdid()`, `mdid()`, `ddid()`,
  `panel_qtt()`, `unc_qte()` (with `target = "qtt"` for the QTT; `lou_qtt()`
  for staggered adoption) and `autoplot()` instead.
* Removed the exported helper `panelize.data()`, which only served those
  functions.
* Removed the exported internal `QTEparams()` parameter container. `unc_qte()`
  now has self-contained estimation and bootstrap code and no longer uses it.
* `formula.tools`, `data.table` and `parallel` are no longer dependencies.
* Requires `ptetools (>= 1.1.0)` and `BMisc (>= 1.4.11)`.
* `unc_qte()`'s covariate argument `xformla` is renamed `xformula`, matching
  the other estimators. `xformla` still works for now, with a deprecation
  warning. The new `clustervars` argument comes right after `weightsname`,
  so calls passing `probs` or later arguments by position need updating.
* `est_method = "or"` (outcome regression) in `unc_qte()` and `lou_qtt()` is
  renamed `"reg"`, the name `did` uses, and the name in the new `cdf_did()`.
  `"or"` still works, with a deprecation warning.

## New features

* `cic()`, `ddid()`, `qdid()`, `mdid()`, `lou_qtt()` and `panel_qtt()` gain a
  `clustervars` argument, as in `did::att_gt()`, to cluster standard errors
  at a level above the unit (ptetools #6). With the empirical bootstrap,
  whole clusters are resampled.
* `unc_qte()` also gains `clustervars`: its bootstrap then resamples whole
  clusters instead of rows.

* New `cdf_did()` (and its workhorse `cdf_did_gt()`): the QTT and ATT under
  parallel trends in the *distribution functions* of untreated potential
  outcomes (Roth and Sant'Anna 2023), with staggered adoption, for panel
  data or repeated cross sections. The counterfactual distribution is
  computed exactly on the data's support and rearranged into a valid
  distribution function; the ATT equals the usual DiD ATT. Unlike the
  rank- and quantile-based estimators, it is identified for discrete
  outcomes. Covariates enter through `est_method = c("aipw", "ipw", "reg")`:
  a doubly robust, an inverse-probability-weighted, and an
  outcome-regression (distribution regression) version of the conditional
  assumption (Sant'Anna and Zhao 2020).

* The estimators now warn, once per call, when the outcome looks discrete:
  more than 10% of observations at one value, or at most 20 distinct values
  (#17). Quantile treatment effects usually require a continuous outcome:
  quantiles jump at mass points, rank-based identification only gives
  bounds, and bootstrap inference can be unreliable. The warning has class
  `qte_discrete_outcome`, the thresholds are the options
  `qte.mass_point_tol` and `qte.discrete_n_unique`, and the diagnostics are
  returned as `outcome_discreteness`. `cic_gt()` no longer warns about
  "extreme ranks" once per group-time cell and bootstrap draw.

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
* `unc_qte()`'s `cband = TRUE` (the default) now reports a uniform
  confidence band, `qte` +/- `c` * `qte.se`, built the same way as the panel
  estimators' uniform bands. Previously `cband` was ignored: `qte.lower` /
  `qte.upper` were always pointwise intervals, though `summary()` labelled
  them "Simult.". The bootstrap critical value `c` was already computed and
  stored. Bands are typically wider under the default; use `cband = FALSE` for
  pointwise intervals. `summary()` and `autoplot()` label and draw the band
  correctly, including for `QTE` objects saved by earlier versions.
* `unc_qte()` with `weightsname`: each bootstrap draw now resamples units
  together with their weights. Previously the weights stayed in their
  original row order while the rows were resampled, so standard errors and
  confidence bands with weights were wrong. Point estimates are unchanged.
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
  columns. Requires `ptetools (>= 1.1.0)`. Results for other column names
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

* With `weightsname`, the panel estimators' overall and event-study results
  (ATT and QTT) now combine cohorts by their weighted shares, as `did` does.
  They used unweighted counts of units, so results were wrong whenever the
  average weight differed across cohorts (ptetools #23; for
  `panel_qtt(pre_copula = "long")`, which aggregates in qte, #13). Integer
  weights now give the same estimates as duplicating units. Results without
  weights are unchanged.

* `mdid()` with covariates no longer fails with "object 'pre_ctrl' not
  found". The outcome regressions looked up their weights where `lm()` could
  not see them, so `mdid()` errored for any non-empty `xformula`.

## Minor improvements

* Aggregated QTTs are now computed exactly: the group-time cells' CDFs are
  mixed from every support point and the mixture inverted exactly, instead of
  on a 1,000-point grid of the pooled outcome (#19; ptetools #25 for the
  estimators that use ptetools' aggregation, and `panel_qtt()`'s
  long-copula aggregation in qte). QTT estimates change slightly (by a few
  thousandths on the package's tests and `mpdta`), and weighted data now
  gives the same QTTs as the equivalent row-replicated data. Requires
  `BMisc (>= 1.4.11)`.
* `unc_qte()` now checks up front that `dname` is a column of `data` coded
  0/1, with an informative error. A treatment coded e.g. 1/2 previously
  failed inside the propensity-score fit with
  "y values must be 0 <= y <= 1" (#1).
* `unc_qte()` drops rows with missing values in the outcome, treatment,
  weights or covariates (with a warning) instead of failing. Its internals
  were rewritten as standalone code; results without weights are unchanged.
* `cic()` and `qdid()` with covariates and `weightsname` (panel data): the
  ATT paired the counterfactual outcomes, sorted by unit id, with weights in
  data row order (#14). ptetools sorts panel data by id, so results were not
  affected in practice; the ATT now uses the same id-ordered weights as the
  counterfactual distribution.
* Invalid sampling weights in the panel estimators (missing, negative, or
  varying over time within a unit) now give a clear error, from ptetools
  (#24), instead of `NaN` or silently inconsistent estimates.
* Repeated cross sections (#18): the help pages and README now say which
  estimators accept `panel = FALSE` (`cic()`, `qdid()`, `mdid()`,
  `cdf_did()`, and `lou_qtt()` with `lagged_outcome_cov = FALSE`) and why
  `ddid()` and `panel_qtt()` need panel data. New tests check that, on panel
  data, `panel = FALSE` gives the panel estimates, and that the estimators
  work on genuine repeated cross sections.
* New tests check that integer weights match row replication for every
  panel estimator, and that rescaling the weights changes nothing (#15).
* `unc_qte()`'s `pscore.reg` no longer carries copies of the data. The fit
  kept the estimator's working data twice over, through `$data`/`$model` and
  through its formula's environment, so a saved `QTE` object for
  `lalonde.psid` shrinks from about 3.6 MB to 0.4 MB. `coef()`, `vcov()`,
  `summary()`, `fitted()` and `predict(newdata = )` are unchanged;
  `anova()` and `update()`, which refit the model, no longer work on it.

## Performance

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

  Combined with the faster empirical bootstrap in ptetools 1.1.0,
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
