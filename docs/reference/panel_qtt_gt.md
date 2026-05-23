# Panel QTT: group-time estimator (Callaway-Li 2019)

Computes the panel QTT and counterfactual outcome distribution for a
single (g,t) cell using the three-period copula-stability estimator of
Callaway and Li (2019). Serves directly as the `attgt_fun` argument to
[`ptetools::pte`](https://rdrr.io/pkg/ptetools/man/pte.html).

**Identification.** Under copula stability
\$\$C(F\_{Y^{\text{pre2}}\|D=1},\\ F\_{\Delta Y^{\text{pre}}\|D=1}) =
C(F\_{Y^{\text{pre1}}(0)\|D=1},\\ F\_{\Delta
Y^{\text{post}}(0)\|D=1})\$\$ and distributional parallel trends on
changes \$\$F\_{\Delta Y^{\text{post}}(0)\|D=1} = F\_{\Delta
Y^{\text{ctrl}}}\$\$ the counterfactual outcome for each treated unit
\\i\\ is \\kcf_i = L_i + C_i,\\ where \\L_i =
Q\_{Y^{\text{pre1}}\|D=1}(u_i)\\ is the treated pre1 quantile at rank
\\u_i = F\_{Y^{\text{pre2}}\|D=1}(Y^{\text{pre2}}\_i)\\ (the Rosenblatt
transform), and \\C_i = Q\_{\Delta Y^{\text{ctrl}}}(v_i)\\ is the
control change at rank \\v_i = F\_{\Delta Y^{\text{pre}}\|D=1}(\Delta
Y^{\text{pre}}\_i)\\. The goal is distributional: \\\\kcf_i\\\\ is a
sample from \\F\_{Y(0)^{\text{post}}\|D=1}\\, not individual
counterfactuals.

**Panel data required for both groups.** Three periods for treated
(pre2, pre1, post) and two for control (pre1, post).

## Usage

``` r
panel_qtt_gt(gt_data, xformula = ~1, ...)
```

## Arguments

- gt_data:

  A `gt_data_frame` with columns `name` (`"pre2"`, `"pre1"`, `"post"`),
  `D`, `Y`, `id`, `.w`, and any covariate columns.

- xformula:

  One-sided formula for covariates. Default `~1`.

- ...:

  Additional arguments passed through by `ptetools`.

## Value

A
[`ptetools::attgt_noif`](https://rdrr.io/pkg/ptetools/man/attgt_noif.html)
object with `att` and, in `extra_gt_returns`, `F0` (counterfactual
ECDF), `F1` (observed post ECDF), and `Fte` (individual-effect ECDF).

## References

Callaway, Brantly and Tong Li. “Quantile Treatment Effects in
Difference-in-Differences Models with Panel Data.” Quantitative
Economics 10(4), pp. 1579-1618, 2019.
