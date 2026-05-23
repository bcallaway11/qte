# Mean Difference-in-Differences: group-time estimator

Computes the MDiD ATT and counterfactual outcome distribution for a
single 2x2 (pre/post x treated/control) data subset. Serves directly as
the `attgt_fun` argument to
[`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html). Panel vs.
repeated cross sections is detected automatically from whether the same
unit ids appear in both periods.

**Identification.** MDiD assumes that the counterfactual distribution of
untreated potential outcomes for the treated group in the post period is
a location shift of the treated group's pre-period distribution:
\$\$F\_{Y(0),\text{post}\|D=1}(y) = F\_{Y,\text{pre}\|D=1}(y -
\Delta)\$\$ where \\\Delta = E\[Y\_{\text{post}}\|D=0\] -
E\[Y\_{\text{pre}}\|D=0\]\\ is the mean DiD (the change in mean outcomes
for the untreated group). The counterfactual quantile function is
therefore \\Q\_{Y(0),\text{post}\|D=1}(\tau) =
Q\_{Y,\text{pre}\|D=1}(\tau) + \Delta\\.

**Covariate adjustment.** With covariates, the scalar shift \\\Delta\\
is replaced by a unit-specific conditional mean shift \\\Delta(X_i) =
E\[Y\_{\text{post}}\|D=0, X_i\] - E\[Y\_{\text{pre}}\|D=0, X_i\]\\,
estimated by weighted OLS on the control group in each period. The
counterfactual for treated pre-period unit \\i\\ is
\\Y\_{\text{pre},i} + \Delta(X_i)\\. The unconditional counterfactual
distribution is the empirical CDF of these shifted values, which by the
law of iterated expectations consistently estimates
\\F\_{Y(0),\text{post}\|D=1}(y) = \int F\_{Y,\text{pre}\|D=1,X=x}(y -
\Delta(x))\\dF\_{X\|D=1}(x)\\.

## Usage

``` r
mdid_gt(gt_data, xformula = ~1, ...)
```

## Arguments

- gt_data:

  A data frame (typically a `gt_data_frame` from `ptetools`) with
  columns `name` (`"pre"` or `"post"`), `D` (treatment dummy), `Y`
  (outcome), `id` (unit identifier), `.w` (sampling weights), and any
  covariate columns referenced by `xformula`.

- xformula:

  One-sided formula for covariates. Default `~1` uses no covariates.
  With covariates, separate weighted OLS regressions are fit on the
  control group in each period; see Details above.

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

Thuysbaert, Bram. “Distributional Comparisons in Difference in
Differences Models.” Working Paper, 2007.
