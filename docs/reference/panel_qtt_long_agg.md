# QTT Aggregation for pre_copula = "long"

Custom `aggregation_fun` for
[`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html) that
handles the two-step weight renormalization needed when
`pre_copula = "long"`: some (g,t) cells are structurally uncomputable
(their pre2 does not exist in the data), and those cells must be
excluded from aggregation with proper weight redistribution.

**Two-step weights (overall):**

1.  Identify surviving groups (those with at least one valid
    post-treatment cell); renormalize \\p_g\\ among those groups only.

2.  Within each surviving group, assign equal weight \\p_g /
    \|\text{valid cells for } g\|\\ to each valid post-treatment cell.

**Dynamic (event-study) weights:** for each event time e, renormalize
\\p_g\\ among groups where (g, g+e) is valid.

**Error behaviour:** if a structurally valid cell returns `NA` att, the
function errors immediately — such NAs indicate a bug in `panel_qtt_gt`,
not a structural limitation.

## Usage

``` r
panel_qtt_long_agg(attgt.list, ptep, extra_gt_returns)
```

## Arguments

- attgt.list:

  List of per-cell attgt results from `compute.pte`.

- ptep:

  ptetools parameter environment.

- extra_gt_returns:

  List of per-cell extra returns (F0, F1, Fte ECDFs).

## Value

A list with `overall_results`, `dyn_results`, `group_results`,
`attgt_results`, `F0_overall`, `F1_overall` — matching the structure
returned by
[`ptetools::qtt_pte_aggregations`](https://rdrr.io/pkg/ptetools/man/qtt_pte_aggregations.html).
