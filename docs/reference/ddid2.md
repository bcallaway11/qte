# ddid2

**Deprecated.** Use
[`ddid`](https://bcallaway11.github.io/qte/reference/ddid.md) instead.

`ddid2` is a legacy two-period wrapper for the Distributional
Difference-in-Differences estimator. The modern replacement
[`ddid`](https://bcallaway11.github.io/qte/reference/ddid.md) uses a
`yname`/`gname`/`tname` interface and supports staggered treatment
adoption natively.

## Usage

``` r
ddid2(
  formla,
  xformla = NULL,
  t,
  tmin1,
  tname,
  data,
  panel = TRUE,
  dropalwaystreated = TRUE,
  idname = NULL,
  probs = seq(0.05, 0.95, 0.05),
  iters = 100,
  alp = 0.05,
  method = "logit",
  se = TRUE,
  retEachIter = FALSE,
  seedvec = NULL,
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

  logical; whether data are panel (must be `TRUE`).

- dropalwaystreated:

  logical; whether to drop always-treated units.

- idname:

  name of the unit id column.

- probs:

  quantile grid.

- iters:

  bootstrap iterations.

- alp:

  significance level.

- method:

  propensity score estimation method (`"logit"`).

- se:

  logical; whether to compute bootstrap SEs.

- retEachIter:

  logical; return per-iteration results.

- seedvec:

  optional vector of seeds for reproducible bootstrap.

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

[`ddid`](https://bcallaway11.github.io/qte/reference/ddid.md)
