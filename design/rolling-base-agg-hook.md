# Design: `aggregation_fun` Hook in ptetools for Rolling-Base Aggregation

**Status:** Planned — not yet implemented  
**Date drafted:** 2026-05-21  
**Affects:** `ptetools` (change), `qte/R/panel_qtt.R` (new aggregation function)

---

## Motivation

The `panel_qtt()` estimator implements the Callaway-Li (2019) panel QTT
identification strategy. The baseline variant ("fixed base") always uses the
two periods immediately preceding treatment onset as the copula-transfer base:
pre2 = g-2, pre1 = g-1 regardless of the event-time horizon e in ATT(g, g+e).

A second variant ("rolling base") uses pre2 = g-e-2, pre1 = g-1 for ATT(g, g+e),
so that the copula is transferred over the same horizon as the treatment effect
being estimated. This is more demanding: it requires e+2 pre-periods, so the
set of computable cells depends on e:

  (g, g+e) is computable  iff  g - e - 2 >= min(tlist)

Consequently, for rolling base:
- Groups with large g can compute effects at all horizons.
- Groups with small g can only compute effects at small e.
- Some (g,t) cells that ATT aggregation would include are simply not estimable.

The current ptetools aggregation path (`attgt_pte_aggregations`,
`qtt_pte_aggregations`, `qott_pte_aggregations`) is selected by a hard-coded
`gt_type` dispatch inside `pte_emp_boot()` (empirical_bootstrap.R:57-63 and
120-126). There is no hook for a caller-supplied aggregation function. Adding
one enables `panel_qtt()` to supply its own aggregation logic — specifically,
two-step weight renormalization — without modifying the core ptetools
aggregation functions (which would risk breaking other estimators).

---

## Proposed Change to ptetools

### 1. Add `aggregation_fun` argument to `pte_emp_boot()`

In `ptetools/R/empirical_bootstrap.R`, change the dispatch block (lines ~56-63
and ~120-126) from:

```r
if (gt_type == "qtt") {
  aggte <- qtt_pte_aggregations(attgt.list, ptep, extra_gt_returns)
} else if (gt_type == "qott") {
  aggte <- qott_pte_aggregations(attgt.list, ptep, extra_gt_returns)
} else {
  aggte <- attgt_pte_aggregations(attgt.list, ptep)
}
```

to:

```r
if (!is.null(aggregation_fun)) {
  aggte <- aggregation_fun(attgt.list, ptep, extra_gt_returns)
} else {
  # Deprecated path: aggregation_fun will be required in a future version.
  warning("aggregation_fun not specified; defaulting based on gt_type. ",
          "This default will be removed in a future version.",
          call. = FALSE)
  if (gt_type == "qtt") {
    aggte <- qtt_pte_aggregations(attgt.list, ptep, extra_gt_returns)
  } else if (gt_type == "qott") {
    aggte <- qott_pte_aggregations(attgt.list, ptep, extra_gt_returns)
  } else {
    aggte <- attgt_pte_aggregations(attgt.list, ptep)
  }
}
```

The same change applies in both the point-estimate block and the bootstrap
loop so that bootstrap draws also use the custom aggregation. `aggregation_fun`
defaults to `NULL` (with a deprecation warning) so existing callers continue
to work during the transition.

### 2. Thread `aggregation_fun` through `pte()`

`pte()` in `ptetools/R/pte.R` is the user-facing entry point. It must accept
`aggregation_fun = NULL` and forward it to `pte_emp_boot()`. No other change
is needed there.

### Expected signature for `aggregation_fun`

```r
aggregation_fun(attgt.list, ptep, extra_gt_returns)
```

This matches the signatures of `qtt_pte_aggregations` and
`qott_pte_aggregations`, so a caller could also pass one of those functions
directly if desired. The return value must match what those functions return
(a named list with at least `overall_results` and `attgt_results`, the same
structure that `pte_emp_boot()` currently gets from the built-in functions).

---

## New Function in qte: `panel_qtt_long_agg()`

This function implements the two-step weight renormalization for `pre_copula = "long"`
aggregation and lives in `qte/R/panel_qtt.R`.

### Two-step weight structure

**Overall ATT:**

1. Determine which groups have at least one computable (g,t) cell (i.e., at
   least one e where g-e-2 >= min(tlist)). Restrict pg to those groups and
   renormalize.
2. Within each surviving group, assign weight pg_g / |valid_cells_for_g| to
   each computable (g,t) cell. Groups with more valid cells do not get more
   total weight than groups with fewer.

**Event-study at horizon e:**

For each e, identify groups where (g, g+e) is computable (g-e-2 >= min(tlist))
and (g+e) is in the observed time periods. Renormalize pg among those groups
for that e. This is equivalent to what `did`'s `balance_e` does at the group
level — entire groups are either in or out for a given e — so no within-group
split is needed.

### Determining validity upfront

The validity of (g, g+e) under `pre_copula = "long"` depends only on g, e,
and min(tlist), which are all known before the bootstrap:

  valid(g, t)  iff  g − (t − g) − 2 ≥ min(tlist)

`panel_qtt_long_agg` computes this mask from `ptep` and uses it in two ways:
- **Invalid cells** (pre2 out of range): excluded from aggregation; weights
  renormalized over the valid set via the two-step procedure.
- **Valid cells that return NA**: hard error — an NA here is unexpected and
  indicates a bug, not a structural limitation.

Do NOT use `is.na(att)` as a proxy for invalidity: that would silently swallow
bugs in `panel_qtt_gt` that produce spurious NAs in structurally valid cells.

### Sketch

```r
panel_qtt_long_agg <- function(attgt.list, ptep, extra_gt_returns) {
  # extract tlist, glist, pg from ptep
  # compute validity mask from structural criterion:
  #   valid(g, t) iff g - (t - g) - 2 >= min(tlist)
  # error if any valid cell has NA att (unexpected, indicates a bug)
  # step 1: surviving groups = those with >= 1 valid cell; renormalize pg
  # step 2: within each surviving group, weight = pg_g / |valid cells for g|
  # call combine_ecdfs restricted to valid cells with renormalized weights
  # return same structure as qtt_pte_aggregations
}
```

The exact internals mirror `qtt_pte_aggregations` (lines 352-530 of
empirical_bootstrap.R) but replace its weight vector with the two-step
renormalized version.

---

## Backward Compatibility

- `aggregation_fun = NULL` triggers a deprecation warning and falls back to
  the current `gt_type`-based dispatch; existing callers continue to work.
- `qtt_pte_aggregations`, `attgt_pte_aggregations`, `qott_pte_aggregations`
  are untouched.
- `panel_qtt(..., pre_copula = "long")` (the new default) passes
  `aggregation_fun = panel_qtt_long_agg` internally.
- `panel_qtt(..., pre_copula = "short")` passes `aggregation_fun =`
  the standard `qtt_pte_aggregations`, making the choice explicit even for
  the simpler case.

---

## Files to Touch

| File | Change |
|------|--------|
| `ptetools/R/empirical_bootstrap.R` | Add `aggregation_fun = NULL` to `pte_emp_boot()`; deprecation warning + dispatch through it in both point-estimate block and bootstrap loop |
| `ptetools/R/pte.R` | Add `aggregation_fun = NULL` to `pte()`; forward to `pte_emp_boot()` |
| `qte/R/panel_qtt.R` | Add `panel_qtt_long_agg()`; add `pre_copula = "long"` argument to `panel_qtt()`; pass `aggregation_fun` accordingly |

---

## Open Questions / Decisions

1. ~~**Argument name:** `base = c("fixed", "rolling")` or `delta_type`?~~
   **Decided (2026-05-21):** `pre_copula = c("long", "short")`. "long" extends
   the copula window to match the event horizon (pre2 = g−e−2); "short" always
   anchors at treatment onset (pre2 = g−2). Default is `"long"`.
2. ~~**combine_ecdfs robustness:** Separately, `combine_ecdfs` in BMisc should
   validate that every list element is callable (error or warning on NULL).~~
   **Done (2026-05-21):** Added `inherits(x, "ecdf")` check at the top of
   `combine_ecdfs`; errors with a message identifying the bad index and its
   class. `make_dist` (and therefore `weighted_ecdf`) already tags return
   values with class `"ecdf"`, so the check is semantically correct.
3. ~~**pre-test cells:**~~ **Resolved (2026-05-21):** Pre-test cells (tp < g)
   always use the sliding base (pre2 = tp−2, pre1 = tp−1) regardless of
   `pre_copula`, since there is no event horizon e to match. These cells are
   guaranteed computable by `required_pre_periods = 2` and land in the
   "expected valid" bucket under both `pre_copula` values — no special casing
   needed in the aggregation.
