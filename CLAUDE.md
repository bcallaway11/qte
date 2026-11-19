# CLAUDE.md — qte Package

## Project Overview

`qte` is a CRAN R package for computing Quantile Treatment Effects (QTE) and
Quantile Treatment Effects on the Treated (QTT). Version 2.0.0 (2026-07) was
a major rewrite that added a unified staggered-treatment-adoption API (built
on `ptetools`). Version 2.1.0 (in development, on `master`) removed the
deprecated two-period wrappers and added distribution regression
(`cond_dist = "dr"`) plus several result-changing bug fixes (see NEWS.md).

**Version:** 2.1.0 (dev); 2.0.0 on CRAN. **GitHub:** https://github.com/bcallaway11/qte.
**pkgdown site:** generated in `docs/` from `_pkgdown.yml` — do not add
top-level `*.md` files casually; pkgdown auto-publishes any top-level `.md`
file except `README.md`/`NEWS.md`/`LICENSE.md`/a few GitHub templates, with
no config override. `CLAUDE.md`/`AGENTS.md` are excluded from `docs/` via a
manual post-build cleanup step (see Development Workflow below).
**CRAN status:** qte 2.0.0 and `ptetools` 1.0.1 are on CRAN. The 2.1.0
release checklist is `dev/TODO.md`.

All estimators use `did`/`ptetools`-style named args (`yname`, `tname`,
`idname`, `biters`, `weightsname`, `cl`, ...):
- `unc_qte` — cross-sectional QTE/QTT under unconfoundedness (`dname` = 0/1
  treatment). Standalone code, no ptetools (below); the staggered analogue
  is `lou_qtt(lagged_outcome_cov = FALSE)`.
- Six panel estimators on `ptetools::pte()` (canonical order): `cic`, `qdid`,
  `panel_qtt`, `ddid`, `mdid`, `lou_qtt` (lagged-outcome unconfoundedness;
  the others are DiD-type), with `gname`, `gt_type`, staggered adoption.

Argument-name wrinkle: `unc_qte` takes `xformla`, the panel estimators take
`xformula` (open item in `dev/TODO.md`).

## Key Architecture

### Current pipeline — ptetools-based (the six panel estimators)
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

### Conditional distributions with covariates (`cond_dist`)
All covariate-adjusted estimators take `cond_dist = c("qr", "dr")` and go
through a dispatch layer in `R/distreg-helpers.R` (`cond_fit`, `cond_cdf`,
`cond_quantile`, `cond_quantile_matrix`): `"qr"` = quantile regression
(helpers in `R/quantile-helpers.R`: rearranged, linearly interpolated
conditional quantiles/CDFs), `"dr"` = distribution regression (logit of
1{Y <= y} on X over a grid). QR results are unchanged by the dispatch layer.
Design notes and open questions: `dev/distreg_report.md`.

### `unc_qte` — standalone cross-sectional estimator
```
unc_qte() → unc_qte_data() → unc_qte_est() (full sample)
          → unc_qte_boot() (resample rows, unc_qte_est() per draw) → QTE
```
- Deliberately independent of ptetools: it has no group/time structure, and
  its `target = "qte"` (population) parameter doesn't fit ptetools'
  effect-on-the-treated aggregation/output classes. Decision 2026-11-19.
- `unc_qte_data()` builds the working data (`.y`, `.d`, `.w` + covariates,
  complete cases); `unc_qte_est()` does IPW / OR / AIPW for QTE or QTT;
  `unc_qte_boot()` is a plain empirical bootstrap via `pbapply::pblapply()`.
- **`QTE`** (qte-class.R) — S3 class for `unc_qte` output, with
  `summary`/`print`/`plot`/`autoplot` methods.
- The legacy `QTEparams`/`setupData`/`bootstrap`/`computeSE`/`SE` machinery
  was deleted in v2.1.0 when `unc_qte` was rewritten (identical estimates; it also fixed the bootstrap,
  which had kept weights in their original row order while resampling rows).

### Key S3 classes
- `QTE` — `unc_qte` output (constructor still exported).
- `pte_results`, `pte_qtt`, `pte_emp_boot` — defined in `ptetools`; returned
  by the six panel estimators.

## File Map

| File | Content |
|------|---------|
| `R/unc_qte.R` | `unc_qte` and its internals `unc_qte_data`, `unc_qte_est`, `unc_qte_boot` |
| `R/cic.R` | `cic_gt`, `cic` (panel or repeated cross-sections via `panel`) |
| `R/qdid.R` | `qdid_gt`, `qdid` |
| `R/panel_qtt.R` | `panel_qtt_gt`, `panel_qtt`, `panel_qtt_long_agg` (custom `aggte_fun` for `pre_copula = "long"`), `three_period_subset` |
| `R/ddid.R` | `ddid_gt`, `ddid` (covariates: controls reweighted by p(X)/(1 − p(X))) |
| `R/mdid.R` | `mdid_gt`, `mdid` |
| `R/lou_qtt.R` | `lou_gt`, `lou_qtt` |
| `R/distreg-helpers.R` | distribution regression (`dr_*`) and the `cond_*` QR/DR dispatch layer |
| `R/quantile-helpers.R` | `wquant`, `wcdf_at`, cached/rearranged QR helpers (`cached_rq`, `qr_cdf`, `qr_quantile`, ...) |
| `R/qte-class.R` | `QTE` class, `summary`/`print`/`plot.QTE` |
| `R/helper-functions.R` | `autoplot.QTE`, `keep_used_columns` |
| `R/data.R` | bundled datasets: `lalonde`, `lalonde.exp`, `lalonde.exp.panel`, `lalonde.psid`, `lalonde.psid.panel` |
| `R/imports.R` | package-level `@import` tags, `globalVariables(".w")` |

`two_by_two_rcs_subset` (repeated cross-sections) lives in `ptetools`, not here.

## Known Issues & Debt

The working list for the next release is `dev/TODO.md`; design history and
rationale are in `dev/NOTES.md`. Highlights:

- `QTE()` (the constructor) is exported though users shouldn't need it.
- `unc_qte` (`xformla`) vs panel estimators (`xformula`) argument naming.
- Parallel bootstrap: a numeric `cl` runs serially on Windows — root cause is
  in ptetools (bcallaway11/ptetools#11).

## Development Workflow

```r
devtools::load_all()    # reload after edits
devtools::document()    # rebuild man/ and NAMESPACE
devtools::test()        # run testthat (~120 expectations; 3 slow bootstrap
                         # tests skipped unless R_SLOW_TESTS=true)
lintr::lint("R/foo.R")  # lint before committing (lint_file() is gone in current lintr)
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

### Vignette caching — do not skip this when editing vignettes

All three vignettes (`intro.qmd`, `panel-estimators.qmd`,
`staggered-adoption.qmd`) precompute their estimator results and cache them
in `vignettes/precomputed/*.rds`. Each vignette's setup chunk does:
```r
cache_file <- "precomputed/<name>-results.rds"
use_cache <- file.exists(cache_file)
if (use_cache) cached <- readRDS(cache_file)
```
and each estimator chunk is `if (use_cache) <- cached$x else { compute; }`,
with a final `save-cache` chunk that writes the whole list once if the cache
didn't already exist. This exists because CRAN rejected the first 2.0.0
submission for **tarball size >10MB** and **vignette rebuild time >10 min**
(see `dev/CRAN_SUBMISSION_PLAN.md`) — without caching, every `R CMD check`
re-runs the full `biters=100` bootstrap for every estimator in every
vignette.

**Two different "cache" concepts here — don't confuse them:**
- `vignettes/precomputed/*.rds` — the deliberate cache above. Tracked in
  git, **not** excluded by `.Rbuildignore`/`.gitignore` — it must ship in
  the tarball so CRAN's build loads it instead of recomputing.
- `vignettes/*_cache/`, `vignettes/*_files/` — quarto/knitr's own automatic
  rendering cache/assets. These *are* excluded via `.Rbuildignore` and
  `.gitignore` (added 2026-07-21) and must never be committed — committing
  `intro_cache/` (9.6MB) is exactly what caused the tarball-size rejection.

**When you edit a vignette's R code** (new estimator call, changed
arguments/data, anything that changes what gets computed), its cache is now
stale:
1. Delete `vignettes/precomputed/<vignette>-results.rds`.
2. Install the package first — `devtools::install()`, not just
   `load_all()`. The quarto vignette engine shells out to a fresh R
   subprocess that only sees the normally-installed package, not an
   in-session `load_all()`.
3. Render from within `vignettes/` so the relative cache path resolves:
   `cd vignettes && quarto render <file>.qmd`. This leaves a stray
   `<file>.html` (and sometimes `<file>.rmarkdown`) next to the source —
   delete them; both are also in `.Rbuildignore`/`.gitignore` since `R CMD
   build` includes anything in `vignettes/` not explicitly excluded there,
   regardless of git status.
4. The full `biters=100` bootstrap can be slow to regenerate. To speed it
   up, temporarily add `cl = 10` (or similar) to the estimator calls, then
   remove it before committing — every estimator function uses `cl` for
   parallel bootstrap iterations.

**Gotcha:** `ptetools` is a hard dependency but isn't loaded by `library(qte)`
alone — only when code actually calls into it (e.g. `cic()`). Loading a
cached result skips that call, so `autoplot()`/`summary()` on a
`pte_qtt`/`pte_results` object can silently dispatch to the wrong method.
`staggered-adoption.qmd`'s setup chunk has `library(ptetools)` to guard
against this; add it to any vignette that plots or summarizes those classes.

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
