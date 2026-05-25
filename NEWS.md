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

* `lou_qte()` — Lagged-outcome unconfoundedness QTT estimator for staggered
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

* Estimators now pass an explicit `aggregation_fun` to `ptetools::pte()`,
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
