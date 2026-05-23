# CLAUDE.md — qte Package

## Project Overview

`qte` is a CRAN R package for computing Quantile Treatment Effects (QTE) and
Quantile Treatment Effects on the Treated (QTET). It covers three main
identification strategies:

1. **Random assignment** — `ci.qte`, `ci.qtet` (Firpo 2007)
2. **Selection on observables** — `DR` (doubly robust)
3. **DiD-based** — `CiC` (Athey-Imbens 2006), `QDiD`, `MDiD`, `ddid2`
   (Callaway-Li 2019, Callaway-Li-Oka 2018), `panel.qtet`

Multi-period / staggered treatment:  `cic2` + `cic_attgt` (beta, built on
`ptetools`).

**CRAN version:** 1.4.0  
**GitHub:** https://github.com/bcallaway11/qte  
**pkgdown site:** generated in `docs/`

## Key Architecture

### Central data-flow pattern
All two-period estimators share this pipeline:
```
User call (e.g., CiC) → QTEparams → setupData → compute.* → bootstrap → QTE
```

- **`QTEparams`** (QTEparams.R) — collects all user arguments into a named list
- **`setupData`** (setupData.R) — extracts treated/untreated subsets and ECDFs
  and **injects them into the caller's environment** via `assign(…, envir=env)`
  where `env <- parent.frame()`. This side-effectful pattern is intentional
  but fragile; handle it carefully.
- **`compute.*`** functions — the statistical workhorse for each method; called
  both for the point estimate and inside each bootstrap iteration
- **`bootstrap`** / **`bootiter`** (bootstrap.R) — empirical bootstrap of SE
- **`QTE`** (qte.R) — S3 class constructor holding estimates and CIs

### Multi-period (staggered) path
```
cic2 → ptetools::pte2 → ptetools::two_by_two_subset → cic_attgt → attgt_noif
```
`cic_attgt` uses `tidyr::pivot_wider` internally. Both `ptetools` and `tidyr`
must be listed in DESCRIPTION (currently missing).

### Key S3 classes
- `QTE` — output of all two-period estimators
- `QTEparams` — parameter container
- `BoundsObj` — output of `bounds`

## Files & Responsibilities

| File | Content |
|------|---------|
| `R/CiC.R` | `compute.CiC`, `CiC`, `cic_attgt`, `cic2` |
| `R/qte.R` | `QTE` class, `panel.qtet`, `compute.panel.qtet` |
| `R/QTEparams.R` | `QTEparams` constructor |
| `R/setupData.R` | `setupData` — injects vars into caller env |
| `R/bootstrap.R` | `bootstrap`, `bootiter` |
| `R/ddid2.R` | `ddid2`, `compute.ddid2` |
| `R/QDiD.R` | `QDiD`, `compute.QDiD` |
| `R/MDiD.R` | `MDiD`, `compute.MDiD` (does NOT use setupData) |
| `R/ciqte.R` | `ci.qte`, `compute.ci.qte` |
| `R/ciqtet.R` | `ci.qtet`, `compute.ci.qtet` |
| `R/DR.R` | `DR` |
| `R/spatt.R` | `spatt`, `compute.spatt` |
| `R/bounds.R` | `bounds`, `BoundsObj` |
| `R/helper-functions.R` | `qtes2mat`, `ggqte`, `diffQ` |
| `R/imports.R` | package-level `@import` tags |
| `R/panel-checks.R` | `panel.checks`, `panelize.data` |
| `R/zzz.R` | `.onLoad` / startup |

## Known Issues & Debt

### Dependency issues (critical)
- `ptetools` and `tidyr` are used in `CiC.R` but absent from DESCRIPTION
- `pte` has been **renamed** to `ptetools`; legacy `pte` emits a deprecation
  warning. All references should move to `ptetools`.

### Documentation
- `cic2`'s `ret_quantile` param has an unfinished inline note:
  "...would be better to return the overall distribution..."
- Many exported internal (`compute.*`, `setupData`) lack `@keywords internal`

### Tests
- **No test suite exists.** `tests/` directory is absent entirely.

### Code quality
- None of the R files have the standard file header block (required per
  global CLAUDE.md conventions)
- `compute.MDiD` does its own data-setup code instead of using `setupData`
- Several internal functions are exported unnecessarily in NAMESPACE
- `setupData` uses `assign(..., envir=parent.frame())` — side-effectful;
  treat carefully during refactors

## Staggered Treatment Adoption (ptetools integration)

The `cic2` function is the multi-period/staggered-treatment version of CiC.
It is in **beta**. The goal is to eventually support:
- `qdid2` — staggered QDiD (not yet implemented)
- `ddid22` — staggered ddid2 (not yet implemented)
- `panel.qtet2` — staggered panel QTET (not yet implemented)

The `attgt` function pattern (`cic_attgt`) is the hook for `ptetools::pte2`.
See `ptetools::did_attgt` for reference signature.

## Development Workflow

```r
devtools::load_all()    # reload after edits
devtools::document()    # rebuild man/ and NAMESPACE
devtools::test()        # run testthat (once tests exist)
lintr::lint_file("R/foo.R")  # lint before committing
rcmdcheck::rcmdcheck()  # full CRAN check
```

Always run `devtools::document()` after editing roxygen2 tags.

## Related Packages

- **ptetools** — multi-period/staggered treatment framework this package
  hooks into (replaces the old `pte` package)
- **BMisc** — utility functions (Brant's own package), heavily used
- **did** — Callaway-Sant'Anna DiD; `mpdta` dataset useful for testing
- **quantreg** — quantile regression (used for covariate adjustment)
