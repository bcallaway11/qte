# Panel QTT (Callaway-Li 2019)

Computes Quantile Treatment Effects on the Treated (QTT) and the Average
Treatment Effect on the Treated (ATT) using the three-period panel
identification strategy of Callaway and Li (2019). Handles two-period
(one post + two pre-treatment periods) and staggered treatment adoption
via `ptetools`. Requires panel data for both groups.

**Identification.** Under a copula stability assumption and
distributional parallel trends on changes, the counterfactual
distribution \\F\_{Y(0)^{\text{post}}\|D=1}\\ is recovered from three
periods of panel data without requiring rank invariance.

**pre_copula.** Controls which pre-treatment periods are used as the
copula-transfer base for post-treatment cells (g, g+e):

- `"long"` (default):

  pre2 = 2\*(g - anticipation) - (g+e) - 2, matching the copula window
  to the event horizon e. Cells where pre2 does not exist in the data
  are excluded from aggregation with two-step weight renormalization.

- `"short"`:

  pre2 = g - anticipation - 2, always anchored at treatment onset.

Pre-test placebo cells always use the sliding base regardless of
`pre_copula`.

## Usage

``` r
panel_qtt(
  yname,
  gname,
  tname,
  idname,
  data,
  xformula = ~1,
  weightsname = NULL,
  control_group = "notyettreated",
  anticipation = 0,
  alp = 0.05,
  cband = TRUE,
  biters = 100,
  cl = 1,
  gt_type = "att",
  probs = NULL,
  pre_copula = "long"
)
```

## Arguments

- yname:

  Name of the outcome variable.

- gname:

  Name of the treatment group variable (first treatment period; 0 for
  never-treated units).

- tname:

  Name of the time period variable.

- idname:

  Name of the unit id variable (required).

- data:

  A data frame.

- xformula:

  One-sided formula for covariates. Default `~1`.

- weightsname:

  Name of the sampling weights column. Default `NULL`.

- control_group:

  `"notyettreated"` (default) or `"nevertreated"`.

- anticipation:

  Number of periods of anticipation. Default `0`.

- alp:

  Significance level. Default `0.05`.

- cband:

  Logical; uniform confidence band if `TRUE` (default).

- biters:

  Number of bootstrap iterations. Default `100`.

- cl:

  Number of parallel clusters. Default `1`.

- gt_type:

  `"att"` (default) or `"qtt"`.

- probs:

  Quantile grid for `gt_type = "qtt"`. Default `seq(0.05, 0.95, 0.05)`.

- pre_copula:

  `"long"` (default) or `"short"`. See Details.

## Value

For `gt_type = "att"`, a `pte_emp_boot` object. For `gt_type = "qtt"`, a
`pte_qtt` object.

## References

Callaway, Brantly and Tong Li. “Quantile Treatment Effects in
Difference-in-Differences Models with Panel Data.” Quantitative
Economics 10(4), pp. 1579-1618, 2019.

## See also

[`cic`](https://bcallaway11.github.io/qte/reference/cic.md),
[`qdid`](https://bcallaway11.github.io/qte/reference/qdid.md),
[`mdid`](https://bcallaway11.github.io/qte/reference/mdid.md),
[`ddid`](https://bcallaway11.github.io/qte/reference/ddid.md)

## Examples

``` r
# \donttest{
data(mpdta, package = "did")

## Panel QTT with rolling pre-period copula (default pre_copula = "long")
res <- panel_qtt(yname = "lemp", gname = "first.treat", tname = "year",
                 idname = "countyreal", data = mpdta,
                 gt_type = "qtt", probs = seq(0.1, 0.9, 0.1), biters = 20)
summary(res)
#> 
#> Overall QTT Curve:
#>  Quantile     QTT Std. Error 95% CB Lower 95% CB Upper
#>       0.1 -0.0270     0.1148      -0.4780       0.4239
#>       0.2 -0.0852     0.0813      -0.4042       0.2339
#>       0.3 -0.0250     0.0630      -0.2722       0.2223
#>       0.4 -0.0680     0.0676      -0.3335       0.1976
#>       0.5 -0.0965     0.0562      -0.3172       0.1242
#>       0.6  0.0070     0.0574      -0.2186       0.2325
#>       0.7  0.0409     0.0386      -0.1106       0.1924
#>       0.8 -0.0201     0.0379      -0.1687       0.1286
#>       0.9 -0.0832     0.1037      -0.4902       0.3238
#> 
# }
```
