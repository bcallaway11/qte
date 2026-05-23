# Quantile Difference-in-Differences: group-time estimator

Computes the QDiD ATT and counterfactual outcome distribution for a
single 2x2 (pre/post x treated/control) data subset. Serves directly as
the `attgt_fun` argument to
[`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html). Panel vs.
repeated cross sections is detected automatically from whether the same
unit ids appear in both periods.

## Usage

``` r
qdid_gt(gt_data, xformula = ~1, ...)
```

## Arguments

- gt_data:

  A data frame (typically a `gt_data_frame` from `ptetools`) with
  columns `name` (`"pre"` or `"post"`), `D` (treatment dummy), `Y`
  (outcome), `id` (unit identifier), `.w` (sampling weights), and any
  covariate columns referenced by `xformula`.

- xformula:

  One-sided formula for covariates. Default `~1` uses no covariates.
  With covariates, conditional quantile regressions are used following
  Athey and Imbens (2006).

- ...:

  Additional arguments passed through by `ptetools`; not used directly.

## Value

A
[`ptetools::attgt_noif`](https://rdrr.io/pkg/ptetools/man/attgt_noif.html)
object with the ATT estimate and, in `extra_gt_returns`, three objects:
`F1` (weighted ECDF of observed treated outcomes in the post period),
`F0` (weighted ECDF of counterfactual untreated outcomes for the treated
group), and `Fte` (weighted ECDF of individual treatment effects; `NULL`
for repeated cross sections).

## References

Athey, Susan and Guido Imbens. “Identification and Inference in
Nonlinear Difference-in-Differences Models.” Econometrica 74(2), pp.
431-497, 2006.
