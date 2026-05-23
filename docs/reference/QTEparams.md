# QTEparams

`QTEparams` is an object that contains all the parameters passed to QTE
methods

## Usage

``` r
QTEparams(
  formla,
  xformla = NULL,
  t = NULL,
  tmin1 = NULL,
  tmin2 = NULL,
  tname = NULL,
  data,
  panel = FALSE,
  w = NULL,
  idname = NULL,
  probs,
  alp = NULL,
  method = NULL,
  plot = NULL,
  se = NULL,
  iters = NULL,
  retEachIter = NULL,
  bootstrapiter = NULL,
  seedvec = NULL,
  pl = NULL,
  cores = NULL,
  est_method = NULL,
  target = NULL
)
```

## Arguments

- formla:

  Should be some y on treatment variable

- xformla:

  a formula for the other covariates such as ~ x1 + x2

- t:

  The last period (not always used)

- tmin1:

  The last pre-treatment period (not always used)

- tmin2:

  The 2nd to last pre-treatment period (not always used)

- tname:

  The name of the column containing time periods (not always used)

- data:

  The name of the data frame

- panel:

  Whether or not the data is panel

- w:

  Additional (usually sampling) weights

- idname:

  The name of the id column used with panel data (not always used)

- probs:

  Which quantiles to produce quantile treatment effects for

- alp:

  The significance level

- method:

  The method to compute the propensity score

- plot:

  boolean for whether or not to plot qtes

- se:

  boolean whether or not to compute standard errors

- iters:

  The number of bootstrap iterations to use to compute standard errors

- retEachIter:

  boolean whether or not to return the full results from each bootstrap
  iteration

- bootstrapiter:

  Used internally for determining whether or not a call is part of
  computing standard errors via the bootstrap

- seedvec:

  A seed to compute the same bootstrap standard errors each time the
  method is called (not always used)

- pl:

  Boolean for whether or not computing bootstrap standard errrors in
  parallel

- cores:

  The number of cores to use if computing standard errors in in parallel

- est_method:

  Estimation method: `"ipw"` (propensity score reweighting), `"or"`
  (outcome regression via quantile regression inversion), or `"aipw"`
  (doubly-robust combination of IPW and OR).

- target:

  Target parameter: `"qte"` (population QTE) or `"qtt"` (QTT, quantile
  treatment effect on the treated).

## Value

A named list with class `"QTEparams"` containing all supplied arguments.
Used internally by legacy estimators (`CiC`, `QDiD`, `MDiD`, `ddid2`,
`panel.qtet`, `ci.qtet`) to pass parameters to `setupData` and
`compute.*` functions.
