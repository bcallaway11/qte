# CDF Difference-in-Differences for a Single 2x2 Cell

Computes the ATT and the distributions needed for the QTT for one
group-time cell under parallel trends in the distribution functions of
untreated potential outcomes. Used as the `attgt_fun` in
[`pte`](https://rdrr.io/pkg/ptetools/man/pte.html) by
[`cdf_did`](https://bcallaway11.github.io/qte/reference/cdf_did.md).

## Usage

``` r
cdf_did_gt(gt_data, xformula = ~1, est_method = c("aipw", "ipw", "reg"), ...)
```

## Arguments

- gt_data:

  A `gt_data_frame` for one (g,t) cell, as produced by
  [`ptetools::two_by_two_subset`](https://rdrr.io/pkg/ptetools/man/two_by_two_subset.html)
  (panel) or
  [`ptetools::two_by_two_rcs_subset`](https://rdrr.io/pkg/ptetools/man/two_by_two_rcs_subset.html)
  (repeated cross sections), with columns `Y`, `D`, `name`, `id`, `.w`
  and any covariates.

- xformula:

  One-sided formula of covariates. Default `~1`.

- est_method:

  With covariates: `"aipw"` (default, doubly robust), `"ipw"`, or
  `"reg"` (outcome regression). Ignored without covariates.

- ...:

  Additional arguments passed by `ptetools`; ignored.

## Value

A
[`ptetools::attgt_noif`](https://rdrr.io/pkg/ptetools/man/attgt_noif.html)
object with the ATT and, in `extra_gt_returns`, `F0` (rearranged
counterfactual distribution), `F1` (treated post-period distribution),
`Fte = NULL`, and `cdf_violation`.

## Details

The counterfactual distribution of untreated outcomes for the treated
group in the post period is \$\$F\_{0}(y) = F\_{\text{pre},D=1}(y) +
F\_{\text{post},D=0}(y) - F\_{\text{pre},D=0}(y),\$\$ evaluated exactly
at every point of the pooled support of the three distributions on the
right. The raw \\F_0\\ need not be a valid distribution function (it can
decrease, or leave \[0, 1\]); it is made one by monotone rearrangement
(sorting its values; Chernozhukov, Fernandez-Val and Galichon 2010) and
clipping to \[0, 1\]. The largest change this makes is returned as
`cdf_violation`; large values are evidence against the identifying
assumption.

The ATT is the mean of the treated group's observed post-period outcomes
minus the mean of the *raw* \\F_0\\, which is what the assumption
identifies. Without covariates this equals the usual DiD ATT.

With covariates, the comparison group's change in CDFs is averaged over
the treated group's covariate distribution. For each \\y\\ this is a DiD
for the outcome \\1\\Y \le y\\\\, estimated as in Sant'Anna and Zhao
(2020): by outcome regression (distribution regression of \\1\\Y \le
y\\\\ on the covariates, in each period), by inverse probability
weighting, or doubly robust (AIPW, the default). IPW is exact on the
data's support points; outcome regression and AIPW use the distribution
regression's threshold grid (see
[`cic_gt`](https://bcallaway11.github.io/qte/reference/cic_gt.md)'s
`cond_dist = "dr"`). The propensity score is a weighted logit of
treatment on the pre-period covariates.

The distribution of individual treatment effects is not identified, so
`Fte` is `NULL`.

## References

Roth, Jonathan, and Pedro H. C. Sant'Anna. “When Is Parallel Trends
Sensitive to Functional Form?” Econometrica 91(2), pp. 737-747, 2023.

Sant'Anna, Pedro H. C., and Jun Zhao. “Doubly Robust
Difference-in-Differences Estimators.” Journal of Econometrics 219(1),
pp. 101-122, 2020.

Chernozhukov, Victor, Ivan Fernandez-Val, and Alfred Galichon. “Quantile
and Probability Curves Without Crossing.” Econometrica 78(3), pp.
1093-1125, 2010.
