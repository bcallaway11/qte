# AGENTS.md — qte Package

> Codex counterpart to CLAUDE.md. Shared conventions live here;
> tool-specific differences are noted where they exist.

## Project Overview

`qte` is a CRAN R package (v2.0.0) for computing Quantile Treatment Effects
(QTE) and Quantile Treatment Effects on the Treated (QTT). Version 2.0.0
(2026-07) was a major rewrite adding a unified staggered-treatment-adoption
API (built on `ptetools`); most of the original two-period estimators are
now deprecated wrappers.

GitHub: https://github.com/bcallaway11/qte. pkgdown site generated in
`docs/` — pkgdown auto-publishes any top-level `.md` file except
`README.md`/`NEWS.md`/`LICENSE.md`/a few GitHub templates, with no config
override, so `CLAUDE.md`/`AGENTS.md` get swept in on every rebuild. Delete
`docs/{AGENTS,CLAUDE}.{md,html}` after each `pkgdown::build_site()` call
(see Development Commands).

CRAN status: archived 2026-07-16; `ptetools` 1.0.1 (hard prerequisite)
resubmitted 2026-07-20, qte resubmission follows once accepted. See
`dev/CRAN_SUBMISSION_PLAN.md`.

Two calling conventions coexist:
1. **Current API** (staggered-capable, `did`/`ptetools`-style: `yname`,
   `gname`, `tname`, `idname`, `biters`, `gt_type`, ...): `unc_qte`
   (cross-sectional, unconfoundedness) plus six DiD-based estimators —
   `cic`, `qdid`, `mdid`, `ddid`, `panel_qtt`, `lou_qtt`.
2. **Legacy formula-based API** (`formla`, `xformla`, `t`, `tmin1`, `iters`,
   `se`, ...): `CiC`, `QDiD`, `MDiD`, `ddid2`, `panel.qtet`, `ci.qte`,
   `ci.qtet` — retained only as `.Deprecated()` wrappers. Planned removal
   ~2026-08 (see dev/NOTES.md).

## Key Architecture

### Current pipeline — ptetools-based
```
User call (e.g. cic()) → attgt_fun (e.g. cic_gt) + subset_fun + aggte_fun
  → ptetools::pte() → pte_results (ATT) / pte_qtt (QTT curve)
```
- Each estimator supplies its own per-group-time-cell function as
  `attgt_fun` (`cic_gt`, `qdid_gt`, `mdid_gt`, `ddid_gt`, `panel_qtt_gt`,
  `lou_gt`), a `subset_fun` (`two_by_two_subset` for panel,
  `two_by_two_rcs_subset` for repeated cross-sections), and an `aggte_fun`
  picked by `gt_type` (`"att"` → `attgt_pte_aggregations`; `"qtt"` →
  `qtt_pte_aggregations`, or `panel_qtt_long_agg` for
  `panel_qtt(pre_copula = "long")`).
- `aggte_fun` was renamed from `aggregation_fun` on 2026-07-19 to avoid a
  parameter-name collision with `contdid`. Must be passed explicitly by
  every caller — `pte()`'s `...` swallows a mismatched name silently rather
  than erroring, so naming drift across `qte`/`ptetools` fails silently.
- `ptetools::pte()` owns the group-time loop, empirical bootstrap, and
  aggregation dispatch. Returns `pte_results` or `pte_qtt` (S3 classes from
  `ptetools`, with their own print/summary/autoplot methods).

### Legacy pipeline — formula-based (deprecated wrappers only)
```
User call → QTEparams → setupData → compute.* → bootstrap → QTE
```
- `setupData` injects extracted subsets/ECDFs into the caller's environment
  via `assign(..., envir=parent.frame())`. Fragile by design; only the
  legacy wrappers still use it. Do not refactor without understanding this
  side effect.
- `compute.*` functions run for both the point estimate and each bootstrap
  draw.
- `compute.MDiD` does its own data setup instead of calling `setupData` —
  known inconsistency, not a bug.

### S3 Classes
- `QTE`, `QTEparams` — legacy output/params (deprecated path only).
- `pte_results`, `pte_qtt`, `pte_emp_boot` — from `ptetools`; returned by
  all six current DiD estimators.

## File Map

| File | Content |
|------|---------|
| `R/unc_qte.R` | `unc_qte`, `compute.unc_qte`; deprecated `ci.qte`, `ci.qtet` |
| `R/cic.R` | `cic_gt`, `cic`; deprecated `CiC`, `compute.CiC` |
| `R/qdid.R` | `qdid_gt`, `qdid`; deprecated `QDiD`, `compute.QDiD` |
| `R/mdid.R` | `mdid_gt`, `mdid`; deprecated `MDiD`, `compute.MDiD` (own data setup) |
| `R/ddid.R` | `ddid_gt`, `ddid`; deprecated `ddid2`, `compute.ddid2` |
| `R/panel_qtt.R` | `panel_qtt_gt`, `panel_qtt`, `panel_qtt_long_agg`, `three_period_subset` |
| `R/lou_qtt.R` | `lou_gt`, `lou_qtt` — new estimator, no legacy equivalent |
| `R/qte-class.R` | `QTE` class, `SE`, `summary`/`print`/`plot.QTE`; deprecated `panel.qtet`, `compute.panel.qtet`, `computeSE` |
| `R/QTEparams.R` | `QTEparams` constructor (legacy) |
| `R/setupData.R` | data-setup utility (injects into parent env, legacy) |
| `R/bootstrap.R` | `bootstrap`, `bootiter` (legacy) |
| `R/panel-checks.R` | `panel.checks`, `panelize.data` (legacy, used by `panel.qtet`) |
| `R/helper-functions.R` | `autoplot.QTE`, `ggqte` (deprecated plotting path) |
| `R/data.R` | datasets: `lalonde`, `lalonde.exp`, `lalonde.exp.panel`, `lalonde.psid`, `lalonde.psid.panel` |

## Known Issues & Debt

Full, actively-maintained list in `dev/NOTES.md`. Highlights:

- Legacy infrastructure deletion planned ~2026-08-01: once the seven
  deprecated wrappers are removed, `R/bootstrap.R`, `R/setupData.R`,
  `R/panel-checks.R`, `R/QTEparams.R`, and most of `R/qte-class.R` become
  dead code.
- `formula.tools` used only in `compute.panel.qtet`'s `pscore` branch —
  removed alongside `panel.qtet` rather than replaced now.
- Several internal `compute.*` functions still exported unnecessarily;
  `globalVariables()` in `R/imports.R` is trimmable once legacy wrappers go.

## Development Commands

```r
devtools::load_all()
devtools::document()
devtools::test()             # 83 fast tests; slow bootstrap tests need
                              # R_SLOW_TESTS=true
lintr::lint_file("R/foo.R")
rcmdcheck::rcmdcheck(args = "--as-cran")
```

`ptetools` and `qte` live in sibling repos (`~/Dropbox/ptetools`,
`~/Dropbox/qte`); load both from source when testing changes across that
boundary, since the installed library version may lag:
```r
devtools::load_all("/home/bmc43193/Dropbox/ptetools")
devtools::load_all(".")
```

After `pkgdown::build_site()`:
```r
file.remove(file.path("docs", c("AGENTS.md", "AGENTS.html",
                                  "CLAUDE.md", "CLAUDE.html")))
```

## Related Packages

- **ptetools** — staggered treatment framework; all six DiD estimators hook
  into `pte()`.
- **BMisc** — utility functions (Brant's package).
- **did** — Callaway-Sant'Anna DiD; `mpdta` useful for testing staggered
  adoption.
- **quantreg** — quantile regression for covariate adjustment.
- **csabounds**, **contdid** — downstream/sibling packages also depending
  on `qte`/`ptetools`; check with `revdeplite` before interface-changing
  releases.
