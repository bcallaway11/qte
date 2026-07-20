# Changelog

## qte 2.0.0

This is a major update that extends all estimators to support staggered
treatment adoption. The API has changed substantially — please read this
section before upgrading.

### New functions

- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  — cross-sectional QTE and QTT under unconfoundedness. Supports IPW,
  outcome regression, and doubly robust (AIPW) estimation. Uses a
  `yname`/`dname` interface consistent with the `did`/`ptetools`
  convention. Replaces
  [`ci.qte()`](https://bcallaway11.github.io/qte/reference/ci.qte.md)
  and
  [`ci.qtet()`](https://bcallaway11.github.io/qte/reference/ci.qtet.md),
  which are retained as deprecated wrappers.

- [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) /
  [`cic_gt()`](https://bcallaway11.github.io/qte/reference/cic_gt.md) —
  Change in Changes estimator for staggered treatment adoption (Athey
  and Imbens 2006). Supports panel and repeated cross sections data, and
  returns either ATT (`gt_type = "att"`) or a full QTT curve
  (`gt_type = "qtt"`). Replaces
  [`CiC()`](https://bcallaway11.github.io/qte/reference/CiC-deprecated.md),
  which is retained as a deprecated wrapper.

- [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md) /
  [`qdid_gt()`](https://bcallaway11.github.io/qte/reference/qdid_gt.md)
  — Quantile Difference in Differences estimator for staggered treatment
  adoption (Callaway, Li, and Oka 2018). Replaces
  [`QDiD()`](https://bcallaway11.github.io/qte/reference/QDiD-deprecated.md),
  which is retained as a deprecated wrapper.

- [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) /
  [`mdid_gt()`](https://bcallaway11.github.io/qte/reference/mdid_gt.md)
  — Mean Difference in Differences estimator for staggered treatment
  adoption (Callaway and Li 2019). Replaces
  [`MDiD()`](https://bcallaway11.github.io/qte/reference/MDiD-deprecated.md),
  which is retained as a deprecated wrapper.

- [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) /
  [`ddid_gt()`](https://bcallaway11.github.io/qte/reference/ddid_gt.md)
  — Distributional Difference in Differences estimator for staggered
  treatment adoption (Callaway and Li 2019). Replaces
  [`ddid2()`](https://bcallaway11.github.io/qte/reference/ddid2.md),
  which is retained as a deprecated wrapper.

- [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  /
  [`panel_qtt_gt()`](https://bcallaway11.github.io/qte/reference/panel_qtt_gt.md)
  — Panel QTT estimator using copula transfer across pre-treatment
  periods (Callaway and Li 2019), with staggered treatment adoption
  support. The `pre_copula` argument controls whether a rolling
  (`"long"`, default) or fixed (`"short"`) base period is used to
  transfer the copula. Replaces
  [`panel.qtet()`](https://bcallaway11.github.io/qte/reference/panel.qtet.md),
  which is retained as a deprecated wrapper.

- `lou_qte()` — Lagged-outcome unconfoundedness QTT estimator for
  staggered treatment adoption. Uses the lagged outcome as a control
  variable to identify the QTT. New; no prior equivalent in the package.

### Deprecated functions

The following functions emit a runtime deprecation warning via
[`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) and will be
removed in v2.1 (target: 2026-08). Each points to its modern replacement
in the warning message.

| Deprecated                                                                  | Replacement                                                               |
|-----------------------------------------------------------------------------|---------------------------------------------------------------------------|
| [`ci.qte()`](https://bcallaway11.github.io/qte/reference/ci.qte.md)         | `unc_qte(target = "qte")`                                                 |
| [`ci.qtet()`](https://bcallaway11.github.io/qte/reference/ci.qtet.md)       | `unc_qte(target = "qtt")`                                                 |
| [`CiC()`](https://bcallaway11.github.io/qte/reference/CiC-deprecated.md)    | [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md)             |
| [`QDiD()`](https://bcallaway11.github.io/qte/reference/QDiD-deprecated.md)  | [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md)           |
| [`MDiD()`](https://bcallaway11.github.io/qte/reference/MDiD-deprecated.md)  | [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md)           |
| [`ddid2()`](https://bcallaway11.github.io/qte/reference/ddid2.md)           | [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md)           |
| [`panel.qtet()`](https://bcallaway11.github.io/qte/reference/panel.qtet.md) | [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md) |

### Removed functions

- `computeDiffSE()` — was dead code; never called anywhere in the
  package. Design notes for a future QTT-vs-ATT test are in
  `dev/NOTES.md`.

- `DR()`, `spatt()`, `bounds()` — removed in earlier 1.5.x development
  work.

- `qtes2mat()`, `diffQ()` — formatting helpers tied to the legacy `QTE`
  class with no internal callers.

### Infrastructure

- Full `testthat` test suite added (83 fast tests; 10 slow bootstrap
  tests skipped by default, enabled via `R_SLOW_TESTS=true`). The
  package previously had no tests.

- All new estimators are backed by
  [`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html). The
  `attgt_fun` hook pattern makes it straightforward to add new
  identification strategies without duplicating bootstrap or aggregation
  code.

- Estimators now pass an explicit `aggte_fun` to
  [`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html),
  replacing the earlier implicit dispatch based on `gt_type`.

- `ptetools` added to `Imports` (was missing from `DESCRIPTION` despite
  being used).

- `Hmisc` and `texreg` removed from `Imports` (were unused after earlier
  cleanup).

- License updated from GPL-2 to GPL-3.

------------------------------------------------------------------------

## qte 1.4.0

- Added `cic2()` for Change in Changes with multiple periods and
  variation in treatment timing, using `ptetools` as a backend. This
  function was experimental; it has been superseded by
  [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) in
  v2.0.0.

## qte 1.3.3

- Bug fix: critical values for uniform confidence bands now correctly
  reflect the `alp` argument rather than being stuck at `alp = 0.05`.

## qte 1.3.2

- Significance level is now passed through in the results.
- Added ability to change significance level in
  [`ggqte()`](https://bcallaway11.github.io/qte/reference/ggqte.md) for
  pointwise confidence intervals without re-estimating.

## qte 1.3.1

CRAN release: 2022-09-01

- Updated maintainer contact information.
- Improvements and updates to vignettes.

## qte 1.3.0

CRAN release: 2019-06-10

- Added first-step quantile regression support for
  [`panel.qtet()`](https://bcallaway11.github.io/qte/reference/panel.qtet.md),
  [`QDiD()`](https://bcallaway11.github.io/qte/reference/QDiD-deprecated.md),
  and
  [`CiC()`](https://bcallaway11.github.io/qte/reference/CiC-deprecated.md).
- Added checks that input data is in the correct format, particularly
  for panel data.

## qte 1.2.3

- Now using empirical quantiles in the CiC method.
- Added code for computing uniform confidence bands.
- General improvements to the codebase, including combining shared code
  across methods.

## qte 1.2.2

CRAN release: 2018-04-30

- Bug fixes to the CiC method.

## qte 1.2

- Added
  [`ddid2()`](https://bcallaway11.github.io/qte/reference/ddid2.md) to
  compute quantile treatment effects under a Difference in Differences
  assumption following Callaway and Li (2016).
- Covariates for all methods now specified using formulas (`xformla`
  argument) rather than passing variable names directly.
- Improvements to bootstrapping standard errors.
- Some extensions to computing the ATT via `spatt()`.
- Preliminary version of multi-period ATT moved to the `did` package.
