# CiC

**Deprecated.** Use
[`cic`](https://bcallaway11.github.io/qte/reference/cic.md) instead.

`CiC` is a legacy two-period wrapper for the Change in Changes
estimator. It accepts a formula and explicit `t`/`tmin1` time period
arguments. The modern replacement
[`cic`](https://bcallaway11.github.io/qte/reference/cic.md) uses a
`yname`/`gname`/`tname` interface and supports staggered treatment
adoption natively via `ptetools`.

## Usage

``` r
CiC(
  formla,
  xformla = NULL,
  t,
  tmin1,
  tname,
  data,
  panel = FALSE,
  se = TRUE,
  idname = NULL,
  alp = 0.05,
  probs = seq(0.05, 0.95, 0.05),
  iters = 100,
  pl = FALSE,
  cores = 2,
  retEachIter = FALSE
)
```

## Arguments

- formla:

  formula of the form `y ~ treat`.

- xformla:

  optional one-sided covariate formula.

- t:

  post-treatment period.

- tmin1:

  pre-treatment period.

- tname:

  name of the time column.

- data:

  data.frame.

- panel:

  logical; whether data are panel.

- se:

  logical; whether to compute bootstrap SEs.

- idname:

  name of the unit id column (required when `panel = TRUE`).

- alp:

  significance level.

- probs:

  quantile grid.

- iters:

  bootstrap iterations.

- pl:

  logical; parallel bootstrap.

- cores:

  number of cores.

- retEachIter:

  logical; return per-iteration results.

## Value

`QTE` object.

## See also

[`cic`](https://bcallaway11.github.io/qte/reference/cic.md)
