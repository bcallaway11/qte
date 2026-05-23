# Distributional DiD: group-time estimator

Computes the distributional DiD ATT and counterfactual outcome
distribution for a single 2x2 (pre/post x treated/control) data subset.
Serves directly as the `attgt_fun` argument to
[`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html).

**Identification.** Under distributional parallel trends and a copula
restriction (Callaway, Li, and Oka 2018), the counterfactual outcome for
each control unit \\j\\ is \$\$kcf_j = \Delta Y\_{\text{ctrl},j} +
Q\_{1,\text{pre}}(u_j)\$\$ where \\\Delta Y\_{\text{ctrl},j} =
Y\_{\text{post},j} - Y\_{\text{pre},j}\\ is the observed change for
control unit \\j\\, \\u_j = F\_{0,\text{pre}}(Y\_{\text{pre},j})\\ is
that unit's rank in the control pre-period distribution, and
\\Q\_{1,\text{pre}}\\ is the quantile function of the treated pre-period
distribution. The unconditional counterfactual distribution
\\F\_{Y(0),\text{post}\|D=1}\\ is then the (weighted) empirical CDF of
\\\\kcf_j\\\\.

Unlike CiC, QDiD, and MDiD, the counterfactual is indexed over *control*
units, not treated units. Consequently `F0` and the ATT counterfactual
term are weighted by `w_pre_ctrl`, and no individual treatment effect
distribution (`Fte`) is returned.

**Panel data required.** The estimator needs the actual change \\\Delta
Y\_{\text{ctrl},j}\\ for each control unit, which requires observing the
same units in both periods.

## Usage

``` r
ddid_gt(gt_data, xformula = ~1, ...)
```

## Arguments

- gt_data:

  A data frame (typically a `gt_data_frame` from `ptetools`) with
  columns `name` (`"pre"` or `"post"`), `D` (treatment dummy), `Y`
  (outcome), `id` (unit identifier), `.w` (sampling weights), and any
  covariate columns referenced by `xformula`. Control units must be
  observed in both periods.

- xformula:

  One-sided formula for covariates. Default `~1` uses no covariates.
  With covariates, the unconditional rank \\u_j\\ is replaced by a
  conditional rank estimated via quantile regression on the control
  pre-period (`QR0tmin1`), and the treated pre-period quantile is
  replaced by a conditional quantile (`QR1tmin1`) evaluated at that rank
  and the control unit's own covariate values.

- ...:

  Additional arguments passed through by `ptetools`; not used directly.

## Value

A
[`ptetools::attgt_noif`](https://rdrr.io/pkg/ptetools/man/attgt_noif.html)
object with the ATT estimate and, in `extra_gt_returns`, `F0` (weighted
ECDF of counterfactual outcomes indexed over control units), `F1`
(weighted ECDF of observed treated post-period outcomes), and
`Fte = NULL` (individual treatment effect distribution is not identified
for this estimator).

## References

Callaway, Brantly, Tong Li, and Tatsushi Oka. “Quantile Treatment
Effects in Difference in Differences Models under Dependence
Restrictions and with Only Two Time Periods.” Journal of Econometrics
206(2), pp. 395-413, 2018.
