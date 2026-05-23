# Three-Period Subset for Panel QTT

Custom `subset_fun` for
[`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html) that pulls
three time periods for each (g,t) cell: pre2, pre1 = g - anticipation -
1, and post = tp. How pre2 is chosen depends on `pre_copula`:

- `"long"` (default):

  pre2 = 2\*(g - anticipation) - tp - 2, so the copula-transfer window
  matches the event horizon e = tp - (g - anticipation). Cells where
  this pre2 does not exist in the data return empty data and
  `panel_qtt_gt` returns `NA`.

- `"short"`:

  pre2 = g - anticipation - 2 (fixed at treatment onset regardless of
  e).

For pre-treatment placebo cells (`tp < g - anticipation`), the base
always slides backwards regardless of `pre_copula`: pre2 = tp - 2, pre1
= tp - 1.
[`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
sets `required_pre_periods = 2` so that the minimum tp in `tlist` is the
third period, guaranteeing tp - 2 is always a valid period.

## Usage

``` r
three_period_subset(
  data,
  g,
  tp,
  control_group = "notyettreated",
  anticipation = 0,
  pre_copula = "long",
  ...
)
```

## Arguments

- data:

  A data frame prepared by
  [`ptetools::setup_pte`](https://rdrr.io/pkg/ptetools/man/setup_pte.html).

- g:

  The current group (first treatment period).

- tp:

  The current time period.

- control_group:

  `"notyettreated"` (default) or `"nevertreated"`.

- anticipation:

  Number of periods of anticipation. Default `0`.

- pre_copula:

  `"long"` (default) or `"short"`. See above.

- ...:

  Additional arguments; silently ignored.

## Value

A list with elements `gt_data`, `n1`, `disidx`.
