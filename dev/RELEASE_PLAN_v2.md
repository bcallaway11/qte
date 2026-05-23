# qte v2.0.0 Release Plan

## Why 2.0?

This update is the largest in the package's history:
- Five new estimators with a completely new calling convention (`unc_qte`,
  `cic`, `qdid`, `mdid`, `ddid`, `panel_qtt`, `lou_qte`)
- All legacy two-period estimators deprecated and slated for removal
- ptetools integration throughout — staggered adoption is now first-class
- NAMESPACE cleaned (removed `Hmisc`, `texreg`, `computeDiffSE`)
- Full test suite added (was zero tests before)

A 2.0 version number accurately signals to users that the API has changed
substantially and they should read the NEWS before upgrading.

---

## Step 1 — DESCRIPTION cleanup

**Remove from Imports:**
- `Hmisc` — was used only by `compute.ci.qtet` (now deleted)
- `texreg` — no longer called anywhere in the package

**Keep for now:**
- `formula.tools` — still used inside `compute.panel.qtet` (deprecated pscore
  path). Remove in the 2026-08 legacy cleanup pass along with `panel.qtet`.

**Verify:**
- `tidyr` — NOT in qte's DESCRIPTION and should not be; qte never calls tidyr
  directly, ptetools handles it internally. No action needed in qte.

**Update Description field:** the current one-line description mentions only
Firpo (2007), Athey-Imbens (2006), Callaway-Li (2019), Callaway-Li-Oka (2018).
Add a sentence describing the staggered adoption API and `lou_qte`.

**Change License: GPL-2 → GPL-3.**

**Bump Version: 1.5.0 → 2.0.0.**

---

## Step 2 — tidyr removal from ptetools (separate ptetools commit)

`tidyr::pivot_wider` is called at exactly two sites in
`ptetools/R/attgt_functions.R` (lines 35 and 173). The operation is a simple
long-to-wide reshape (two columns: `pre`, `post`) that is trivially expressible
with `data.table::dcast()`, which ptetools already imports:

```r
# Replace:
gt_data_outcomes <- tidyr::pivot_wider(
  gt_data[, c("D", "id", "period", "name", "Y")],
  id_cols = c(id, D), names_from = name, values_from = Y
)

# With:
gt_data_outcomes <- data.table::dcast(
  data.table::as.data.table(gt_data[, c("D", "id", "period", "name", "Y")]),
  id + D ~ name, value.var = "Y"
)
```

After replacing both sites, remove `tidyr` from ptetools DESCRIPTION Imports
and the `@import tidyr` tag from `ptetools/R/imports.R`. Run ptetools tests to
confirm nothing breaks. Commit and push ptetools separately before the qte
v2.0.0 submission.

---

## Step 3 — NEWS.md rewrite

The current 1.5.0 section is a stub. Rewrite it as 2.0.0 with full coverage:

**New functions:**
- `unc_qte()` — cross-sectional QTE/QTT under unconfoundedness (IPW/OR/AIPW);
  replaces `ci.qte()` and `ci.qtet()`
- `cic()`, `cic_gt()` — CiC for staggered adoption; replaces `CiC()`
- `qdid()`, `qdid_gt()` — QDiD for staggered adoption; replaces `QDiD()`
- `mdid()`, `mdid_gt()` — MDiD for staggered adoption; replaces `MDiD()`
- `ddid()`, `ddid_gt()` — distributional DiD for staggered adoption; replaces `ddid2()`
- `panel_qtt()`, `panel_qtt_gt()` — panel QTT with copula transfer for
  staggered adoption; replaces `panel.qtet()`
- `lou_qte()` — lagged-outcome unconfoundedness QTT for staggered adoption
  (new; no prior equivalent)

**Deprecated (with runtime warnings):**
- `CiC()`, `QDiD()`, `MDiD()`, `ddid2()`, `panel.qtet()`, `ci.qte()`,
  `ci.qtet()` — all retained with `.Deprecated()` messages pointing to
  their replacements. Planned removal: v2.1 (target 2026-08).

**Removed:**
- `computeDiffSE()` — uncalled; design notes in `dev/NOTES.md`
- `DR()`, `spatt()`, `bounds()` — removed in earlier 1.5.x work

**Infrastructure:**
- Full testthat suite added (83 fast tests, 10 slow tests skipped by default)
- ptetools integration: all new estimators backed by `ptetools::pte()`
- `aggregation_fun` hook: estimators pass explicit aggregation functions
  (no more fallback dispatch in ptetools)

---

## Step 4 — Examples audit

For every exported function, check that:

1. **Primary new functions** (`unc_qte`, `cic`, `qdid`, `mdid`, `ddid`,
   `panel_qtt`, `lou_qte`) each have at least one runnable `@examples` block.
   Currently these may have no examples. Add minimal examples using the
   `lalonde` data (cross-sectional) or `mpdta` from `did` (staggered).
   Examples should run in < 5 seconds — keep `biters` low (≤ 20) or use
   `\donttest{}` for any call that bootstraps.

2. **Deprecated wrappers** (`CiC`, `QDiD`, `MDiD`, `ddid2`, `panel.qtet`,
   `ci.qte`, `ci.qtet`) that emit deprecation warnings during CRAN check:
   wrap existing examples in `\dontrun{}` so CRAN does not run them (CRAN
   checks treat unexpected warnings as a NOTE or WARNING depending on the
   check level). Add a one-liner showing the replacement instead:
   ```r
   # See ?cic for the modern replacement
   ```

3. **Internal helpers** (`compute.panel.qtet`, `computeSE`, `QTEparams`,
   `panelize.data`): confirm they are either `@keywords internal` (no examples
   needed) or have minimal working examples. `compute.panel.qtet` is currently
   exported and has its own Rd page — it should get `@keywords internal`.

4. **`QTE` class constructor**: has examples via `panel.qtet` examples.
   Audit whether a direct `QTE(...)` example makes sense or if it should be
   `@keywords internal`.

---

## Step 5 — README rewrite

**Format:** `README.qmd` (delete `README.Rmd`). Output `format: gfm`.

**Structure:**
```
# qte — Quantile Treatment Effects in R
[badges: CRAN, R-CMD-check]

## Overview
One paragraph: what QTEs are, what the package does, two bullet groups
(cross-sectional estimators / DiD-based estimators).

## Installation

## Quick start — selection on observables (lalonde data)
unc_qte() example; show plot with autoplot() (now the standard; ggqte() is
deprecated).

## Staggered treatment adoption (mpdta data)
cic() and/or panel_qtt() example; show event-study plot via autoplot()
(autoplot.pte_results for ATT, autoplot.pte_qtt for QTT).

## Available estimators
Table: function | method | staggered? | data type
(covers unc_qte, cic, qdid, mdid, ddid, panel_qtt, lou_qte)
Brief note that CiC/QDiD/MDiD/ddid2/panel.qtet/ci.qte/ci.qtet are
deprecated wrappers pointing to the above.

## Documentation and vignettes
Link to pkgdown site, list vignettes.
```

**Data choice:** lalonde for the quick-start (familiar, self-contained),
mpdta for staggered. The `did` package is already in Suggests so mpdta is
available without a new dependency.

---

## Step 6 — Vignettes rewrite

Delete: `ddid2.Rmd`, `panel-qtet.Rmd` (both reference deprecated APIs, too
short to update).

Replace `R-QTEs.Rmd` and add new:

### Vignette 1: `intro.qmd`
**Title:** "Quantile Treatment Effects in R"
Replaces `R-QTEs.Rmd`. Content:
- What are QTEs and why do we care (brief, 1-2 paragraphs)
- Random assignment: `unc_qte(target="qte")` with lalonde.exp
- Selection on observables: `unc_qte(target="qte")` with covariates on lalonde.psid
- QTT under selection on observables: `unc_qte(target="qtt")`
- Reading and plotting output: `summary()` and `autoplot()` / `plot()`.
  Use `autoplot(result)` — this calls `autoplot.QTE`, the new standard.
  Do NOT use `ggqte()` (deprecated).

### Vignette 2: `did-estimators.qmd`
**Title:** "DiD-Based Quantile Treatment Effects"

All six DiD-based estimators (`cic`, `qdid`, `mdid`, `ddid`, `panel_qtt`,
`lou_qte`) are conceptually peers and all appear in both DiD vignettes. The
split is by **purpose**, not by method or data structure.

Vignette 2 is the **conceptual and introductory** vignette. Goal: "what is
this method and why would I use it?"
- Brief overview: DiD-based identification of QTEs, how distributional
  parallel trends and related assumptions extend the standard DiD logic
- One section per estimator: identification assumption, key restriction, when
  it is appropriate vs. alternatives. Order: CiC → QDiD → DDiD → MDiD →
  panel_qtt → lou_qte. MDiD and DDiD sections are intentionally shorter (fewer
  distinguishing features relative to CiC/QDiD).
- Minimal code per estimator — just enough to show the call signature and
  basic output (lalonde.psid.panel for CiC/QDiD/MDiD/DDiD; mpdta for
  panel_qtt and lou_qte since they need richer panel structure). Use
  `biters = 0` or very low `biters` to keep compile time short.
- Closing table: estimator | assumption | panel/RCS | two-period/staggered

### Vignette 3: `staggered-adoption.qmd`
**Title:** "Staggered Treatment Adoption with the qte Package"

Vignette 3 is the **applied workflow** vignette. Goal: "here is how to run
the full analysis with staggered adoption data." Uses `mpdta` throughout.
- Brief setup: what staggered adoption is, the `gname`/`tname`/`idname`
  interface shared by all estimators
- Run all six estimators on mpdta; compare overall ATT/QTT results
- Aggregation: overall, group-specific, event-study — show plots for each
  using `autoplot()`. Plotting dispatch by return class:
    - `gt_type = "att"` → `pte_results` → `autoplot.pte_results` (event-study
      via `autoplot(res, type = "dynamic")`)
    - `gt_type = "qtt"` → `pte_qtt` → `autoplot.pte_qtt` (overall QTT curve
      via `autoplot(res)`, event-study via `autoplot(res, type = "dynamic")`)
- ATT vs. QTT: `gt_type = "att"` vs `gt_type = "qtt"` side by side for one
  estimator (e.g. CiC), interpret the difference
- Brief note on `panel_qtt` pre_copula options and `lou_qte` lagged-outcome
  argument since those have extra tuning choices the others don't

**All vignettes:** `.qmd` format; `VignetteEngine: quarto::html`; shared
bibliography `refs.bib`.

---

## Step 7 — _pkgdown.yml audit

The reference section almost certainly still organizes functions around legacy
names and groupings. Reorganize into:

```yaml
reference:
  - title: "Cross-sectional estimators"
    contents: [unc_qte, ci.qte, ci.qtet]

  - title: "DiD-based estimators — staggered adoption"
    contents: [cic, cic_gt, qdid, qdid_gt, mdid, mdid_gt, ddid, ddid_gt,
               panel_qtt, panel_qtt_gt, lou_qte]

  - title: "Deprecated two-period wrappers"
    contents: [CiC, QDiD, MDiD, ddid2, panel.qtet]

  - title: "QTE class and output"
    contents: [QTE, SE, summary.QTE, plot.QTE, ggqte]

  - title: "Data"
    contents: [lalonde, lalonde.exp, lalonde.exp.panel,
               lalonde.psid, lalonde.psid.panel]

  - title: "Internal / legacy infrastructure"
    contents: [QTEparams, panelize.data, computeSE, compute.panel.qtet]
```

---

## Step 8 — R CMD check pass

Run `rcmdcheck::rcmdcheck(args = "--as-cran")` and resolve all ERRORs,
WARNINGs, and NOTEs before submission. Expected issues to check:

- Undocumented arguments in any new function
- Missing `@return` or `@examples` in exported functions
- Any `\dontrun{}` that should be `\donttest{}`
  (CRAN prefers `\donttest{}` for slow examples; `\dontrun{}` for examples
  that genuinely cannot run, e.g., because they require external state)
- NOTE: "no visible binding for global variable" — may appear in functions
  that use NSE or `data.table` syntax; suppress with `globalVariables()` or
  explicit `NULL` assignments at package level
- Check that all packages in Imports are actually called with `::` or imported
  via `@import`/`@importFrom` (this will catch any remaining dead imports after
  the Hmisc/texreg removal)

**Update `cran-comments.md`** with the platform, R version, check results
summary, and a note explaining the deprecation strategy.

---

## Step 9 — pkgdown rebuild

After all of the above:
```r
pkgdown::build_site()
```
Check that:
- All new functions appear in the Reference section
- Vignettes build without errors
- The deprecated-functions section appears clearly
- Plots render correctly in the vignette articles

---

## Step 10 — Merge to master and tag

1. Ensure branch `major-update-2026` passes R CMD check cleanly.
2. Open a PR or merge directly (your call) to `master`.
3. Tag: `git tag v2.0.0`.
4. Submit to CRAN.

---

## Rough order of work within a session

1. DESCRIPTION (removes Hmisc/texreg, GPL-3, version 2.0.0, description text)
2. ptetools: drop tidyr (separate commit/push)
3. NEWS.md full rewrite
4. Examples audit (read each Rd, fix/add examples)
5. README.qmd
6. Vignette 1: intro.qmd
7. Vignette 2: did-estimators.qmd
8. Vignette 3: staggered-adoption.qmd (expand)
9. Delete ddid2.Rmd, panel-qtet.Rmd, README.Rmd
10. _pkgdown.yml reorganization
11. rcmdcheck pass
12. pkgdown build
13. Merge + tag + CRAN submission

---

## Open questions — RESOLVED

- **`cic2()`**: Already removed in earlier cleanup — not present in any R file
  or NAMESPACE. No action needed.
- **`ggqte()`**: Deprecated in favour of `autoplot.QTE` (new, lives in qte).
  `plot.QTE` now calls `autoplot()`. `autoplot.pte_results` and
  `autoplot.pte_qtt` already exist in ptetools. Done.
- **`compute.panel.qtet`**: Already `@keywords internal` with no `@export` —
  not present in NAMESPACE. No action needed.
