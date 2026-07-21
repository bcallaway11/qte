## Resubmission after archival

This package was archived on CRAN on 2026-07-16 as a consequence of its
dependency 'ptetools' being archived (ptetools was archived for calling a
deprecated 'BMisc' function; see ptetools's own cran-comments.md). ptetools
1.0.1 has since been resubmitted, fixes that issue, and is back on CRAN.
This is the follow-up resubmission of qte, now depending on the corrected
`ptetools (>= 1.0.1)`.

This is also a major update (1.4.0 -> 2.0.0):

* New unified staggered treatment adoption API built on 'ptetools' for
  all difference-in-differences estimators (`cic`, `qdid`, `mdid`, `ddid`,
  `panel_qtt`).
* New estimators: `unc_qte()` (cross-sectional QTE/QTT under
  unconfoundedness) and `lou_qtt()` (lagged-outcome unconfoundedness).
* The seven legacy two-period estimators (`CiC`, `QDiD`, `MDiD`, `ddid2`,
  `panel.qtet`, `ci.qte`, `ci.qtet`) are retained as deprecated wrappers
  around the new API and emit a deprecation message; they are scheduled
  for removal in a future release.

See NEWS.md for the full list of changes.

## Test environments

* Local Ubuntu, R 4.x — all checks passed without issues.
* GitHub Actions:
    - Windows-latest (R release)
    - Windows-latest (R devel)
    - macOS-latest (R release)
    - Ubuntu-latest (R release)
    - Ubuntu-latest (R devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE — new submission, package was
  previously archived (expected; see "Resubmission after archival" above).

## Downstream dependencies

Checked reverse dependencies on CRAN (`sccic`) and on GitHub (`csabounds`);
no issues found.
