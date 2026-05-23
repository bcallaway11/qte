# ci.qte

Deprecated. Use
[`unc_qte`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
instead.

## Usage

``` r
ci.qte(
  formla,
  xformla = NULL,
  x = NULL,
  data,
  w = NULL,
  probs = seq(0.05, 0.95, 0.05),
  se = TRUE,
  iters = 100,
  alp = 0.05,
  method = "logit",
  retEachIter = FALSE,
  printIter = FALSE,
  pl = FALSE,
  cores = 2,
  biters = NULL,
  cl = NULL
)
```

## Arguments

- formla:

  formula of the form `y ~ treatment`.

- xformla:

  one-sided covariate formula. Default `NULL`.

- x:

  ignored (deprecated argument).

- data:

  data.frame.

- w:

  sampling weight vector. Cannot be forwarded automatically; use
  `weightsname` in `unc_qte` instead.

- probs:

  quantile grid.

- se:

  ignored; bootstrap is always run in `unc_qte`.

- iters:

  number of bootstrap iterations.

- alp:

  significance level.

- method:

  propensity score model.

- retEachIter:

  ignored.

- printIter:

  ignored.

- pl:

  logical; whether to run bootstrap in parallel.

- cores:

  number of cores.

- biters:

  alias for `iters`; takes precedence if supplied.

- cl:

  alias for cores; takes precedence if supplied.

## Value

A `QTE` object; same structure as
[`unc_qte`](https://bcallaway11.github.io/qte/reference/unc_qte.md).
