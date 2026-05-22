# qte API and File Organization Cleanup Plan

## Motivation

The package currently mixes older two-period formula-based APIs with newer
`ptetools`-based staggered-adoption APIs. This is workable during a transition,
but it makes the package harder to navigate:

- file names sometimes reflect legacy exported functions rather than the main
  current API;
- deprecated wrappers and new user-facing functions live side-by-side without a
  fully consistent documentation strategy;
- tests sometimes exercise old wrappers where the new API should be primary;
- some internal `compute.*` helpers remain visible in documentation even though
  they are only compatibility machinery.

The goal is to make the new APIs primary while preserving backwards
compatibility through clearly deprecated wrappers.

## Deprecation Direction

The following functions should be treated as legacy compatibility wrappers:

- `ci.qte()` -> `unc_qte()`
- `ci.qtet()` -> `unc_qte(target = "qtt")`
- `CiC()` -> `cic()`
- `QDiD()` -> `qdid()`
- `MDiD()` -> `mdid()`
- `ddid2()` -> `ddid()`

`ci.qtet()` should be deprecated explicitly. The suggested message is:

```r
.Deprecated("unc_qte")
```

with documentation explaining that users should call:

```r
unc_qte(..., target = "qtt")
```

for Quantile Treatment Effects on the Treated under unconfoundedness.

## Proposed File Organization

Do not rename files until an explicit old -> new mapping is approved. The
following layout is the target organization.

| Current file | Proposed file | Main contents |
|---|---|---|
| `R/unc_qte.R` | `R/unc_qte.R` | `unc_qte()`, `compute.unc_qte()`, deprecated `ci.qte()` |
| `R/ciqtet.R` | merge into `R/unc_qte.R` or rename to `R/ci_qtet.R` | deprecated `ci.qtet()`, `compute.ci.qtet()` if still needed |
| `R/CiC.R` | `R/cic.R` | `cic_gt()`, `cic()`, deprecated `CiC()`, `compute.CiC()` |
| `R/QDiD.R` | `R/qdid.R` | `qdid_gt()`, `qdid()`, deprecated `QDiD()`, `compute.QDiD()` |
| `R/MDiD.R` | `R/mdid.R` | `mdid_gt()`, `mdid()`, deprecated `MDiD()`, `compute.MDiD()` |
| `R/ddid2.R` | `R/ddid.R` | `ddid_gt()`, `ddid()`, deprecated `ddid2()`, `compute.ddid2()` |
| `R/panel_qtt.R` | `R/panel_qtt.R` | `panel_qtt_gt()`, `panel_qtt()`, QTT aggregation helpers |
| `R/qte.R` | `R/qte-class.R` | `QTE()`, `SE()`, print/summary/plot methods, `computeSE()` |
| `R/qte.R` | `R/data.R` | Lalonde data documentation |

## Documentation Updates

Recommended documentation rules:

- New lower-case / snake-case APIs should have the primary user-facing help
  pages: `unc_qte`, `cic`, `qdid`, `mdid`, `ddid`, `panel_qtt`.
- Legacy wrappers should have short help pages that say they are deprecated and
  point users to the replacement function.
- Internal `compute.*` helpers should either be undocumented or marked clearly
  as internal. They should not be presented as primary user APIs.
- `ci.qtet()` documentation should be rewritten as a deprecation page pointing
  to `unc_qte(target = "qtt")`.

## Test Updates

Recommended test split:

- Main estimator tests should use the current APIs:
  `unc_qte`, `cic`, `qdid`, `mdid`, `ddid`, `panel_qtt`.
- Legacy wrappers should get small compatibility tests only:
  - they emit the expected deprecation warning;
  - they return the same broad object class as before;
  - no extensive numerical regression tests unless needed for CRAN/backward
    compatibility.
- Avoid using deprecated wrappers inside "fast smoke tests" unless the warning
  is explicitly expected.

## Suggested Order of Work

1. Deprecate `ci.qtet()` and update its documentation.
2. Update tests so `unc_qte(target = "qtt")` is the primary QTT-under-
   unconfoundedness test path.
3. Regenerate documentation with `devtools::document()`.
4. Review `NAMESPACE` after documentation generation.
5. Only after tests and docs are stable, consider file renames using an
   approved old -> new mapping.

## Open Questions

- Should deprecated wrappers remain exported for one CRAN release or longer?
- Should `compute.ci.qtet()` be retained as an internal helper, or should
  `ci.qtet()` delegate directly to `unc_qte(target = "qtt")`?
- Should data documentation be moved out of `R/qte.R` into a dedicated
  `R/data.R` file as part of the same cleanup or kept separate?
