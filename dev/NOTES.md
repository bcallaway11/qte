# Developer Notes — qte Package

Running notes on design decisions, API conventions, known issues, and
outstanding work. Not formal documentation — just the things worth remembering
across sessions. See also `CLAUDE.md` for project architecture and known debt.

---

## API & Interface

### Two Interface Styles

The package currently has two distinct calling conventions that coexist:

**1. Legacy formula-based API** (original qte estimators):

```r
formla, xformla, t, tmin1, tname, data, panel, idname,
probs, iters, alp, se, pl, cores, retEachIter
```

**2. `did` / `ptetools`-style API** (newer staggered-treatment functions):

```r
yname, gname, tname, idname, data, xformla, anticipation,
alp, biters, cband, boot_type, weightsname
```

Legacy APIs should remain stable for backward compatibility. All new
staggered/multi-period functions should follow the `did` / `ptetools` style.

### Specific Mismatches to Be Aware Of

| Topic | Legacy qte | did / ptetools | Recommendation |
|---|---|---|---|
| Outcome/group spec | `formla`, `t`, `tmin1` | `yname`, `gname` | New code: use name-based |
| Covariates | `xformla` | `xformla` (qte/did), `xformula` (some ptetools) | Use `xformla` for qte-facing code |
| Bootstrap iters | `iters` | `biters` | New code: `biters` |
| Bootstrap on/off | `se` | `boot_type` | New code: `boot_type` |
| Weights | `w` (vector) | `weightsname` (column name) | New code: `weightsname` |
| Parallelism | `pl`, `cores` | `cl` (ptetools) | New code: `cl` |
| Confidence bands | `alp` only | `cband`, `alp` | New code: both |

`panel` default varies across legacy functions (`CiC`, `QDiD`, `MDiD`, `spatt`
default `FALSE`; `ddid2` defaults `TRUE`; `panel.qtet` is panel-only). Not a
bug, but can surprise users switching between estimators.

`MDiD()` does not expose `pl` / `cores` while neighboring estimators do — minor
inconsistency, low priority to fix.

---

## Design Decisions

### `panel_qtt`: `pre_copula = c("long", "short")`

The Callaway-Li (2019) panel QTT identification strategy requires two
pre-treatment periods to transfer the copula. There are two natural choices:

- **`"short"` (fixed base):** pre2 = g − 2, pre1 = g − 1 for all horizons e.
  Always anchors at treatment onset. Requires only 2 pre-periods regardless of e.
- **`"long"` (rolling base, default):** pre2 = 2g − tp − 2, pre1 = g − 1 for
  ATT(g, g+e). Extends the copula window to match the event horizon — more
  demanding (requires e + 2 pre-periods), so fewer (g, t) cells are computable
  at longer horizons.

**Validity check** uses `pre2_long %in% original_tp` (set membership), not
`>= min(tlist)`, to handle unequally-spaced panels correctly.

**NA handling:** The aggregation function (`panel_qtt_long_agg`) uses a
structural validity mask computed from first principles. Valid cells that
return `NA` trigger a hard error — that indicates a bug in `panel_qtt_gt`,
not a structural limitation. Do NOT use `is.na(att)` as a proxy for
structural invalidity.

**Pre-test cells** always use a sliding base (pre2 = tp − 2, pre1 = tp − 1)
regardless of `pre_copula`, since there is no event horizon to match.

### `aggregation_fun` Hook in ptetools

Added to `pte()`, `panel_empirical_bootstrap()`, and `qtt_empirical_bootstrap()`
(ptetools commit 32d899d, 2026-05-21). When `NULL`, emits a deprecation warning
and falls back to `gt_type`-based dispatch. Expected signature:

```r
aggregation_fun(attgt.list, ptep, extra_gt_returns)
```

`panel_qtt()` uses this hook to pass `panel_qtt_long_agg` (for
`pre_copula = "long"`) or `ptetools::qtt_pte_aggregations` (for `"short"`).

---

## Outstanding Work & Future Ideas

### Next Update — Delete Legacy Infrastructure (Target: 2026-08-01)

**What to keep:** `ci.qte` (cross-sectional Firpo 2007 estimator, no staggered
extension planned). `ci.qtet` will get a new ptetools-based staggered version
before the legacy version is removed — do not delete it until that is done.

Once the remaining legacy two-period functions (`CiC`, `QDiD`, `MDiD`,
`ddid2`, `panel.qtet`, `ci.qtet`) are removed, the following
files become dead code and should be deleted in the same pass:

- `R/bootstrap.R` — `bootiter()` and `bootstrap()` are only called by legacy
  `compute.*` functions.
- `R/setupData.R` — `setupData()` is only called by legacy `compute.*`
  functions and `panel.checks()`.
- `R/panel-checks.R` — `panel.checks()` and `panelize.data()` are only called
  by legacy functions. Note: `panelize.data` is exported, so its removal needs
  a NEWS.md entry.
- `R/QTEparams.R` — the `QTEparams` parameter container is only used by the
  legacy pipeline.
- `R/qte.R` — holds the `QTE` S3 class and `panel.qtet`; can go once all
  `compute.*` callers are removed.

Deletion order: remove the legacy estimator files first, then these
infrastructure files, then regenerate NAMESPACE with `devtools::document()`.

### High Priority

- **Install updated BMisc**: Tests currently require loading BMisc from source
  (`devtools::load_all('/home/bmc43193/Dropbox/BMisc')`) because the installed
  version still uses the old `dflist` parameter in `combine_ecdfs`. Install the
  updated version once BMisc is published or pinned.

- **`DESCRIPTION` dependency cleanup**: `ptetools` and `tidyr` are used but not
  listed in `DESCRIPTION` Imports. `cic2` will fail for users without them.

- **Pass explicit `aggregation_fun` in existing callers**: `cic()`, `qdid()`,
  `mdid()`, `ddid()` currently rely on the deprecated NULL fallback in
  `pte_emp_boot()`, which generates warnings during tests. Pass the appropriate
  function explicitly to silence those.

### Medium Priority

- **Rolling base for other estimators**: `pre_copula` currently only implemented
  for `panel_qtt`. Could add analogous support to `cic_gt`, `qdid_gt`, etc.
  Design doc was at `design/rolling-base-agg-hook.md` (deleted 2026-05-21,
  implementation complete for `panel_qtt`).

- **`cic2` cleanup**: `cic2()` calls `pte2()`, `setup_pte`, and
  `two_by_two_subset` unqualified; NAMESPACE doesn't import ptetools. Fix
  before treating `cic2` as stable.

### Print / Summary / Plot for Quantile Results — Needs Revisiting

The new ptetools-based estimators return `pte_qtt` (quantile curve) and
`pte_emp_boot` (ATT) objects. Need to audit whether print, summary, and plot
methods for these are adequate, and whether they belong in `ptetools` (as
generic infrastructure) or in `qte` (as package-specific wrappers).

`pte_qtt` has basic print/summary in `ptetools/R/classes.R` (added 2026-05-21)
but it was not designed with the full quantile use case in mind. Questions:

- Is the `summary.pte_qtt` output useful as-is, or does it need richer formatting?
- Should there be an `autoplot.pte_qtt` in ptetools, or a `ggqte`-style function in qte?
- Does `pte_emp_boot` need a quantile-aware print/summary path for when
  `gt_type = "qtt"` results are stored in it?

Revisit alongside the plotting discussion below.

### Plotting — Needs Revisiting

`ggqte()` in `R/helper-functions.R` is the only plotting function remaining.
It was written for the legacy `QTE` class and will need to be rethought when
those functions are removed. Questions to resolve:

- Should there be a `plot`/`autoplot` S3 method for `pte_qtt` and `pte_emp_boot`?
- Does `ggqte` have external users worth preserving, or can it be replaced cleanly?
- What should the default plot look like for the new staggered estimators
  (overall curve, event-study, group-level)?

Defer until the legacy cleanup is further along.

### Low Priority / Future Ideas

- **Uniform parallel interface**: `MDiD` lacks `pl`/`cores`; standardize across
  all estimators when touching those files.
- **`panel` defaults**: Consider whether inconsistent defaults across `CiC`,
  `QDiD`, `MDiD`, `ddid2` should be harmonized (breaking change, low urgency).
- **Staggered extensions**: `qdid2`, `ddid22`, `panel.qtet2` are planned but
  not yet implemented (see CLAUDE.md).
