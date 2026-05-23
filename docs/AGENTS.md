# AGENTS.md — qte Package

> Codex counterpart to CLAUDE.md. Shared conventions live here;
> tool-specific differences are noted where they exist.

## Project Overview

`qte` is a CRAN R package (v1.4.0) for computing Quantile Treatment
Effects (QTE) and Quantile Treatment Effects on the Treated (QTET).

Identification strategies covered: 1. **Random assignment** — `ci.qte`,
`ci.qtet` (Firpo 2007) 2. **Selection on observables** — `DR` (doubly
robust) 3. **DiD-based** — `CiC` (Athey-Imbens 2006), `QDiD`, `MDiD`,
`ddid2`, `panel.qtet`

Multi-period / staggered treatment: `cic2` + `cic_attgt` (beta, built on
`ptetools`).

GitHub: <https://github.com/bcallaway11/qte>

## Key Architecture

### Two-period estimator pipeline

    User call → QTEparams → setupData → compute.* → bootstrap → QTE

- **`setupData`** injects variables into the caller’s environment via
  `assign(..., envir=parent.frame())`. This is intentional but fragile.
  Do not refactor without understanding this side-effect.
- **`compute.*`** functions are called both for point estimates and
  inside bootstrap iterations.

### Multi-period (staggered) path

    cic2 → ptetools::pte2 → two_by_two_subset → cic_attgt → attgt_noif

`ptetools` and `tidyr` must appear in DESCRIPTION (currently missing).

### S3 Classes

- `QTE` — all two-period estimator output
- `QTEparams` — parameter container
- `BoundsObj` — output of `bounds`

## File Map

| File                   | Content                                                    |
|------------------------|------------------------------------------------------------|
| `R/CiC.R`              | `compute.CiC`, `CiC`, `cic_attgt`, `cic2`                  |
| `R/qte.R`              | `QTE` class, `panel.qtet`                                  |
| `R/QTEparams.R`        | `QTEparams` constructor                                    |
| `R/setupData.R`        | data-setup utility (injects into parent env)               |
| `R/bootstrap.R`        | `bootstrap`, `bootiter`                                    |
| `R/ddid2.R`            | `ddid2`, `compute.ddid2`                                   |
| `R/QDiD.R`             | `QDiD`, `compute.QDiD`                                     |
| `R/MDiD.R`             | `MDiD`, `compute.MDiD` (has own data setup — inconsistent) |
| `R/ciqte.R`            | `ci.qte`, `compute.ci.qte`                                 |
| `R/ciqtet.R`           | `ci.qtet`, `compute.ci.qtet`                               |
| `R/DR.R`               | `DR`                                                       |
| `R/spatt.R`            | `spatt`, `compute.spatt`                                   |
| `R/bounds.R`           | `bounds`, `BoundsObj`                                      |
| `R/helper-functions.R` | `qtes2mat`, `ggqte`, `diffQ`                               |

## Known Issues & Debt

- `ptetools` and `tidyr` missing from DESCRIPTION
- `pte` package has been renamed to `ptetools`; deprecation warning
  fires on load
- No test suite (`tests/` directory absent)
- No standard file headers on any R files
- `compute.MDiD` uses its own data setup instead of `setupData`
- Several internal `compute.*` functions exported unnecessarily

## Planned Extensions

- `qdid2` — staggered QDiD (not yet implemented)
- `ddid22` — staggered ddid2 (not yet implemented)
- `panel.qtet2` — staggered panel QTET (not yet implemented)

## Development Commands

``` r
devtools::load_all()
devtools::document()
devtools::test()
lintr::lint_file("R/foo.R")
rcmdcheck::rcmdcheck()
```

## Related Packages

- **ptetools** — staggered treatment framework (hooks for `cic2`)
- **BMisc** — utility functions (Brant’s package)
- **did** — Callaway-Sant’Anna DiD; `mpdta` useful for testing
- **quantreg** — quantile regression for covariate adjustment
