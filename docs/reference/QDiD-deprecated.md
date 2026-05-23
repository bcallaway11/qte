# QDiD

**Deprecated.** Use
[`qdid`](https://bcallaway11.github.io/qte/reference/qdid.md) instead.

`QDiD` is a legacy two-period wrapper for the Quantile
Difference-in-Differences estimator. It accepts a formula and explicit
`t`/`tmin1` time period arguments. The modern replacement
[`qdid`](https://bcallaway11.github.io/qte/reference/qdid.md) uses a
`yname`/`gname`/`tname` interface and supports staggered treatment
adoption natively.

## Usage

``` r
QDiD(
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
  retEachIter = FALSE,
  pl = FALSE,
  cores = NULL,
  biters = NULL,
  cl = NULL
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

  name of the unit id column.

- alp:

  significance level.

- probs:

  quantile grid.

- iters:

  bootstrap iterations.

- retEachIter:

  logical; return per-iteration results.

- pl:

  logical; parallel bootstrap.

- cores:

  number of cores.

- biters:

  alias for `iters`; takes precedence if supplied.

- cl:

  alias for `cores`; takes precedence if supplied.

## Value

`QTE` object.

## See also

[`qdid`](https://bcallaway11.github.io/qte/reference/qdid.md)
