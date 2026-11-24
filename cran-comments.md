## Submission summary

This is a minor release (2.0.0 -> 2.1.0). Main changes (see NEWS.md):

* Removes the two-period functions deprecated in 2.0.0 (`CiC()`, `QDiD()`,
  `MDiD()`, `ddid2()`, `panel.qtet()`, `ci.qte()`, `ci.qtet()`, `ggqte()`),
  as announced in that release, plus the exported internals `QTEparams()` and
  `panelize.data()`.
* New features: distribution regression for the covariate-adjusted
  estimators (`cond_dist = "dr"`), and `clustervars` for clustered standard
  errors in all estimators. `unc_qte()`'s `xformla` argument is renamed
  `xformula`; the old name still works, with a deprecation warning.
* Bug fixes that change results, including covariate-adjusted conditional
  quantiles (#9), `ddid()` with covariates (#10), the AIPW ATT (#11), and
  sampling weights throughout (weighted cohort shares, a weighted bootstrap
  that kept weights in row order in `unc_qte()`). `mdid()` with covariates
  no longer errors.
* Quantile regressions now use the Frisch-Newton interior point method
  (`quantreg::rq(method = "fn")`). With a mass point in the outcome, the
  default Barrodale-Roberts simplex could cycle forever and hang R.
* Weighted quantiles in `cic()` and `panel_qtt()` no longer depend on
  floating-point rounding, which made results differ across platforms.
* Requires 'ptetools' (>= 1.1.0).
* Drops the 'data.table', 'parallel' and 'formula.tools' dependencies.

## Follow-up to 'ptetools' 1.1.0

This release requires 'ptetools' (>= 1.1.0), now on CRAN (same
maintainer). The 'ptetools' reverse-dependency check flagged 6 test
failures in 'qte' 2.0.0, all fixed here:

* 2 errors: 'ptetools' 1.1.0 stores its internal columns under new names,
  which the subset functions in 'qte' 2.0.0 did not read; this release
  reads the new names.
* 4 changed values: tests hardcoding median quantile treatment effects,
  which change because 'ptetools' 1.1.0 corrects how group-time QTTs are
  aggregated. This release updates them to the corrected values.

## Test environments

* Local: Ubuntu 26.04.1 LTS, R 4.5.2
* GitHub Actions (all Status: OK, 0 errors | 0 warnings | 0 notes):
    - windows-latest (R release)
    - windows-latest (R devel)
    - macOS-latest (R release)
    - ubuntu-latest (R release)
    - ubuntu-latest (R devel)

## R CMD check results

`rcmdcheck::rcmdcheck(args = "--as-cran")`:

0 errors | 0 warnings | 0 notes

(Locally, there is one NOTE that HTML validation and math rendering in the
manual were skipped because HTML Tidy and V8 are not installed on the check
machine; this is a property of the local setup, not the package.)

## Downstream dependencies

* CRAN: `sccic` 0.1.1 (lists qte in Suggests): OK.
* GitHub only (not on CRAN): `csabounds` could not be checked because its
  dependency `distreg` is not installable; it also uses `qte::CiC()` in its
  examples, which this release removes after its 2.0.0 deprecation. Its
  maintainer is the qte maintainer, and it is not on CRAN.
