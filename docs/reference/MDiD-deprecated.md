# MDiD

**Deprecated.** Use
[`mdid`](https://bcallaway11.github.io/qte/reference/mdid.md) instead.

`MDiD` is a legacy two-period wrapper for the Mean
Difference-in-Differences estimator. The modern replacement
[`mdid`](https://bcallaway11.github.io/qte/reference/mdid.md) uses a
`yname`/`gname`/`tname` interface and supports staggered treatment
adoption natively.

## Usage

``` r
MDiD(
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
  biters = NULL
)
```

## Arguments

- formla:

  formula of the form `y ~ treat`.

- xformla:

  optional one-sided covariate formula (unused in this estimator;
  accepted for consistency with other wrappers).

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

- retEachIter:

  logical; return per-iteration results.

- biters:

  alias for `iters`; takes precedence if supplied.

## Value

`QTE` object.

## See also

[`mdid`](https://bcallaway11.github.io/qte/reference/mdid.md)
