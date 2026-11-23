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
* Requires 'ptetools' (>= 1.1.0).
* Drops the 'data.table', 'parallel' and 'formula.tools' dependencies.

## Coordinated submission with 'ptetools'

This release is submitted together with 'ptetools' 1.1.0 (same
maintainer), which it requires. 'ptetools' 1.1.0 stores its internal
columns under new names, which breaks the subset functions in the current
CRAN version of 'qte' (2.0.0); this release reads the new names. Please
check it against 'ptetools' 1.1.0.

## Test environments

* Local: Ubuntu 26.04.1 LTS, R 4.5.2
* GitHub Actions: TODO
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
