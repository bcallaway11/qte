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
* Requires 'ptetools' (>= 1.0.3), released just before this submission.
* Drops the 'data.table', 'parallel' and 'formula.tools' dependencies.

## Test environments

* Local: TODO (OS, R version)
* GitHub Actions: TODO
    - windows-latest (R release)
    - windows-latest (R devel)
    - macOS-latest (R release)
    - ubuntu-latest (R release)
    - ubuntu-latest (R devel)

## R CMD check results

`rcmdcheck::rcmdcheck(args = "--as-cran")`: TODO

## Downstream dependencies

TODO: run `.revdeplite/run_revdeplite.R`.

* CRAN: `sccic` lists qte only in Suggests and does not call it; expected no
  problems.
* GitHub only (not on CRAN): `csabounds` uses `qte::CiC()` in its examples,
  which this release removes after its 2.0.0 deprecation. Its maintainer is
  the qte maintainer, so no notice is needed; it is not on CRAN.
