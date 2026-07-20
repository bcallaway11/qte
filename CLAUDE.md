# CLAUDE.md — qte Package

## Project Overview

`qte` is a CRAN R package for computing Quantile Treatment Effects (QTE) and
Quantile Treatment Effects on the Treated (QTT). Version 2.0.0 (2026-07) was
a major rewrite that added a unified staggered-treatment-adoption API (built
on `ptetools`) alongside the original two-period estimators, most of which
are now deprecated wrappers.

**Version:** 2.0.0. **GitHub:** https://github.com/bcallaway11/qte.
**pkgdown site:** generated in `docs/` from `_pkgdown.yml` — do not add
top-level `*.md` files casually; pkgdown auto-publishes any top-level `.md`
file except `README.md`/`NEWS.md`/`LICENSE.md`/a few GitHub templates, with
no config override. `CLAUDE.md`/`AGENTS.md` are excluded from `docs/` via a
manual post-build cleanup step (see Development Workflow below).
**CRAN status:** was archived 2026-07-16; `ptetools` 1.0.1 (a hard
prerequisite) was resubmitted 2026-07-20 and qte resubmission follows once
it's accepted. See `dev/CRAN_SUBMISSION_PLAN.md`.

Two calling conventions coexist:
1. **Current API** (staggered-adoption-capable, `did`/`ptetools`-style named
   args: `yname`, `gname`, `tname`, `idname`, `biters`, `gt_type`, ...):
   `unc_qte` (cross-sectional, unconfoundedness) and six DiD-based estimators
   — `cic`, `qdid`, `mdid`, `ddid`, `panel_qtt`, `lou_qtt`.
2. **Legacy formula-based API** (`formla`, `xformla`, `t`, `tmin1`, `iters`,
   `se`, ...): `CiC`, `QDiD`, `MDiD`, `ddid2`, `panel.qtet`, `ci.qte`,
   `ci.qtet` — all retained only as `.Deprecated()` wrappers around the
   current API. Planned removal ~2026-08 (see dev/NOTES.md).

## Key Architecture

### Current pipeline — ptetools-based (all six DiD estimators + implicitly unc_qte's own bootstrap)
```
User call (e.g. cic()) → attgt_fun (e.g. cic_gt) + subset_fun + aggte_fun
  → ptetools::pte() → pte_results (ATT) / pte_qtt (QTT curve)
```
- Each estimator supplies its own per-group-time-cell function as
  `attgt_fun` (e.g. `cic_gt`, `qdid_gt`, `mdid_gt`, `ddid_gt`, `panel_qtt_gt`,
  `lou_gt`), a `subset_fun` (`ptetools::two_by_two_subset` for panel data,
  `two_by_two_rcs_subset` for repeated cross-sections), and an `aggte_fun`
  chosen based on `gt_type` (`"att"` → `ptetools::attgt_pte_aggregations`;
  `"qtt"` → `ptetools::qtt_pte_aggregations`, or `panel_qtt_long_agg` for
  `panel_qtt(pre_copula = "long")`).
- `aggte_fun` was renamed from `aggregation_fun` on 2026-07-19 (ptetools
  commit) to avoid a parameter-name collision with `contdid`. Must be passed
  explicitly by every caller — `ptetools::pte()`'s `...` silently swallows a
  mismatched name instead of erroring, so a naming drift here fails silently.
- `ptetools::pte()` owns the group-time loop, the empirical bootstrap, and
  dispatches to the chosen `aggte_fun`. Results are `pte_results` (scalar
  ATT + event-study) or `pte_qtt` (full QTT curve) — both S3 classes defined
  in `ptetools`, with their own `print`/`summary`/`autoplot` methods.

### Legacy pipeline — formula-based (deprecated wrappers only)
```
User call (e.g. CiC) → QTEparams → setupData → compute.* → bootstrap → QTE
```
- **`setupData`** (setupData.R) injects extracted treated/untreated subsets
  and ECDFs into the caller's environment via `assign(…, envir=parent.frame())`.
  Side-effectful and fragile; only touched by the legacy wrappers now.
- **`compute.*`** functions (`compute.CiC`, `compute.QDiD`, `compute.MDiD`,
  `compute.ddid2`, `compute.panel.qtet`, `compute.unc_qte`) are the
  statistical workhorse, called for both the point estimate and each
  bootstrap draw.
- **`bootstrap`/`bootiter`** (bootstrap.R) — legacy empirical bootstrap.
- **`QTE`** (qte-class.R) — S3 class holding legacy estimates/CIs; `SE()`
  builds its SE sub-object.
- `compute.MDiD` does its own data setup instead of calling `setupData` —
  known inconsistency, not a bug.

### Key S3 classes
- `QTE`, `QTEparams` — legacy two-period estimator output/params (deprecated
  path only).
- `pte_results`, `pte_qtt`, `pte_emp_boot` — defined in `ptetools`; returned
  by all six current DiD estimators.

## File Map

| File | Content |
|------|---------|
| `R/unc_qte.R` | `unc_qte`, `compute.unc_qte`; deprecated wrappers `ci.qte`, `ci.qtet` |
| `R/cic.R` | `cic_gt`, `cic`; deprecated wrapper `CiC`, `compute.CiC` |
| `R/qdid.R` | `qdid_gt`, `qdid`; deprecated wrapper `QDiD`, `compute.QDiD` |
| `R/mdid.R` | `mdid_gt`, `mdid`; deprecated wrapper `MDiD`, `compute.MDiD` (own data setup, doesn't use `setupData`) |
| `R/ddid.R` | `ddid_gt`, `ddid`; deprecated wrapper `ddid2`, `compute.ddid2` |
| `R/panel_qtt.R` | `panel_qtt_gt`, `panel_qtt`, `panel_qtt_long_agg` (custom `aggte_fun` for `pre_copula = "long"`), `three_period_subset` |
| `R/lou_qtt.R` | `lou_gt`, `lou_qtt` — new estimator, no legacy equivalent |
| `R/qte-class.R` | `QTE` class, `SE`, `summary`/`print`/`plot.QTE`; deprecated wrapper `panel.qtet`, `compute.panel.qtet`, `computeSE` |
| `R/QTEparams.R` | `QTEparams` constructor (legacy) |
| `R/setupData.R` | `setupData` — injects vars into caller env (legacy) |
| `R/bootstrap.R` | `bootstrap`, `bootiter` (legacy empirical bootstrap) |
| `R/panel-checks.R` | `panel.checks`, `panelize.data` (legacy, used by `panel.qtet`) |
| `R/helper-functions.R` | `autoplot.QTE`, `ggqte` (deprecated plotting path) |
| `R/data.R` | bundled datasets: `lalonde`, `lalonde.exp`, `lalonde.exp.panel`, `lalonde.psid`, `lalonde.psid.panel` |
| `R/imports.R` | package-level `@import` tags, `globalVariables()` |

## Known Issues & Debt

See `dev/NOTES.md` for the full, actively-maintained list (interface
mismatches between the two calling conventions, `panel_qtt`'s `pre_copula`
design, plotting/print method gaps for `pte_qtt`/`pte_emp_boot`, planned new
estimators). Highlights:

- **Legacy infrastructure deletion planned ~2026-08-01**: once the seven
  deprecated wrappers are removed, `R/bootstrap.R`, `R/setupData.R`,
  `R/panel-checks.R`, `R/QTEparams.R`, and most of `R/qte-class.R` become
  dead code and should be deleted in the same pass.
- `formula.tools` is used only inside `compute.panel.qtet`'s `pscore` branch
  — scheduled for removal alongside `panel.qtet` itself rather than replaced
  now.
- Several `compute.*`/internal functions are still exported unnecessarily in
  NAMESPACE; most of `globalVariables()` in `R/imports.R` becomes trimmable
  once the legacy wrappers are gone.

## Development Workflow

```r
devtools::load_all()    # reload after edits
devtools::document()    # rebuild man/ and NAMESPACE
devtools::test()        # run testthat (83 fast tests; 10 slow bootstrap
                         # tests skipped unless R_SLOW_TESTS=true)
lintr::lint_file("R/foo.R")  # lint before committing
rcmdcheck::rcmdcheck(args = "--as-cran")  # full CRAN check
```

`ptetools` and `qte` are developed in sibling local repos
(`~/Dropbox/ptetools`, `~/Dropbox/qte`). When testing changes that touch the
`ptetools` boundary, load both from source rather than relying on the
installed library version, which may lag:
```r
devtools::load_all("/home/bmc43193/Dropbox/ptetools")
devtools::load_all(".")
```
After `pkgdown::build_site()`, delete the auto-published internal-notes
pages it regenerates every time:
```r
file.remove(file.path("docs", c("AGENTS.md", "AGENTS.html",
                                  "CLAUDE.md", "CLAUDE.html")))
```

Always run `devtools::document()` after editing roxygen2 tags.

## Related Packages

- **ptetools** — multi-period/staggered treatment framework all six DiD
  estimators hook into via `pte()`.
- **BMisc** — utility functions (Brant's own package), heavily used.
- **did** — Callaway-Sant'Anna DiD; `mpdta` dataset useful for testing
  staggered adoption.
- **quantreg** — quantile regression (used for covariate adjustment).
- **csabounds**, **contdid** — downstream/sibling packages that also depend
  on `qte`/`ptetools`; check with `revdeplite` (`.revdeplite/run_revdeplite.R`)
  before releases that change shared interfaces.
