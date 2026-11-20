# lou_gt

Estimates the ATT and counterfactual distribution \\F_0 =
F\_{Y(0)\|D=1}\\ for a single (g,t) cell under lagged-outcome
unconfoundedness. Called internally by
[`lou_qtt`](https://bcallaway11.github.io/qte/reference/lou_qtt.md) via
[`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html).

## Usage

``` r
lou_gt(
  gt_data,
  xformula = ~1,
  lagged_outcome_cov = TRUE,
  est_method = "ipw",
  cond_dist = c("qr", "dr"),
  ...
)
```

## Arguments

- gt_data:

  A `gt_data_frame` from
  [`ptetools::two_by_two_subset`](https://rdrr.io/pkg/ptetools/man/two_by_two_subset.html)
  or
  [`ptetools::two_by_two_rcs_subset`](https://rdrr.io/pkg/ptetools/man/two_by_two_rcs_subset.html).
  Contains columns `Y`, `D`, `name` ("pre"/"post"), `id`, and `.w`.

- xformula:

  One-sided covariate formula. Default `~1`.

- lagged_outcome_cov:

  Logical; if `TRUE` (default), the pre-period outcome is appended to
  `xformula` as an additional covariate. Requires panel data (unit ids
  present in both periods).

- est_method:

  One of `"ipw"`, `"reg"`, or `"aipw"` (`"or"` is a deprecated alias for
  `"reg"`). Default `"ipw"`.

- cond_dist:

  For `est_method = "reg"` or `"aipw"`, how the conditional distribution
  of the untreated outcome is estimated: `"qr"` (default) for quantile
  regression or `"dr"` for distribution regression. See
  [`cic_gt`](https://bcallaway11.github.io/qte/reference/cic_gt.md).
  Ignored for `"ipw"`.

- ...:

  Additional arguments (ignored).

## Value

A list from
[`ptetools::attgt_noif`](https://rdrr.io/pkg/ptetools/man/attgt_noif.html)
with `attgt` (ATT) and `extra_gt_returns` containing `F0` and `F1`
(weighted ECDFs of counterfactual and observed treated outcomes).
