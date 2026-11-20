# Change in Changes: group-time estimator

Computes the Change in Changes (CiC) ATT and counterfactual outcome
distribution for a single 2x2 (pre/post x treated/control) data subset.
Serves directly as the `attgt_fun` argument to
[`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html). Panel vs.
repeated cross sections is detected automatically from whether the same
unit ids appear in both periods.

## Usage

``` r
cic_gt(gt_data, xformula = ~1, cond_dist = c("qr", "dr"), ...)
```

## Arguments

- gt_data:

  A data frame (typically a `gt_data_frame` from `ptetools`) with
  columns `name` (`"pre"` or `"post"`), `D` (treatment dummy), `Y`
  (outcome), `id` (unit identifier), and any covariate columns
  referenced by `xformula`.

- xformula:

  One-sided formula for covariates. Default `~1` uses no covariates.
  With covariates, the control group's conditional outcome distributions
  are estimated by the method in `cond_dist`, following Athey and Imbens
  (2006).

- cond_dist:

  How the conditional distributions are estimated when `xformula` has
  covariates: `"qr"` (default) for quantile regression at the
  percentiles 0.01, ..., 0.99, or `"dr"` for distribution regression (a
  logit of \\1\\Y \le y\\\\ on the covariates at a grid of thresholds;
  Chernozhukov, Fernandez-Val and Melly 2013). Distribution regression
  handles outcomes with mass points (e.g. many zeros) naturally. Ignored
  without covariates.

- ...:

  Additional arguments passed through by `ptetools`; not used directly.

## Value

A
[`ptetools::attgt_noif`](https://rdrr.io/pkg/ptetools/man/attgt_noif.html)
object with the ATT estimate and, in `extra_gt_returns`, three objects:
`F1` (ECDF of observed treated outcomes in the post period), `F0` (ECDF
of counterfactual untreated outcomes for the treated group), and `Fte`
(ECDF of individual treatment effects under rank invariance; `NULL` for
repeated cross sections).

## References

Athey, Susan and Guido Imbens. “Identification and Inference in
Nonlinear Difference-in-Differences Models.” Econometrica 74(2), pp.
431-497, 2006.

Chernozhukov, Victor, Ivan Fernandez-Val and Blaise Melly. “Inference on
Counterfactual Distributions.” Econometrica 81(6), pp. 2205-2268, 2013.
