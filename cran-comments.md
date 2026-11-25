## Submission summary

This version of `qte` (2.1.0)fixes several bugs and improves the package's functionality. First, the current version of `qte` (2.0.0) has errors that are due to upstream changes in `ptetools`. This version fixes those errors.  The main other changes are (see NEWS.md for a full list):

* Removes the two-period functions deprecated in 2.0.0 (`CiC()`, `QDiD()`,
  `MDiD()`, `ddid2()`, `panel.qtet()`, `ci.qte()`, `ci.qtet()`, `ggqte()`),
  as announced in that release, plus the exported internals `QTEparams()` and
  `panelize.data()`.
* Distribution regression for the covariate-adjusted
  estimators (`cond_dist = "dr"`)
* `clustervars` for clustered standard errors in all estimators.
* `unc_qte()`'s `xformla` argument is renamed
  `xformula`; the old name still works, with a deprecation warning.
* Added support for sampling weights throughout
* Quantile regressions now use the Frisch-Newton interior point method
  (`quantreg::rq(method = "fn")`) in stead of `quantreg::rq(method = "br")` which could cycle forever in some cases with mass points in the outcome.
* Requires 'ptetools' (>= 1.1.0).
* Drops the 'data.table', 'parallel' and 'formula.tools' dependencies.

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

## Downstream dependencies

* `sccic` 0.1.1 (lists qte in Suggests): OK.
