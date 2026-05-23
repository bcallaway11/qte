# panel.qtet

**Deprecated.** Use
[`panel_qtt`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
instead. `panel_qtt` supports staggered treatment adoption and exposes
the `pre_copula` option for controlling how the copula is transferred
across pre-treatment periods.

`panel.qtet` computes the Quantile Treatment Effect on the Treated
(QTET) using the method of Callaway and Li (2015). This method should be
used when the researcher wants to invoke a Difference in Differences
assumption to identify the QTET. Relative to the other Difference in
Differences methods available in the `qte` package, this method's
assumptions are more intuitively similar to the identifying assumptions
used in identifying the Average Treatment Effect on the Treated (ATT).

Additionally, this method can accommodate covariates in a more flexible
way than the other Difference in Differences methods available. In order
to accommodate covariates, the user should specify a vector `x` of
covariate names. The user also may specify a method for estimating the
propensity score. The default is logit.

`panel.qtet` can only be used in some situations, however. The method
requires three periods of panel data where individuals are not treated
until the last period. The data should be formatted as a panel; the
names of columns containing time periods and ids for each cross
sectional unit need to be passed to the method.

## Usage

``` r
panel.qtet(
  formla,
  xformla = NULL,
  t,
  tmin1,
  tmin2,
  tname,
  data,
  idname,
  probs = seq(0.05, 0.95, 0.05),
  iters = 100,
  alp = 0.05,
  method = c("qr", "pscore"),
  se = TRUE,
  retEachIter = FALSE,
  pl = FALSE,
  cores = NULL,
  biters = NULL,
  cl = NULL
)
```

## Arguments

- formla:

  The formula y ~ d where y is the outcome and d is the treatment
  indicator (d should be binary), d should be equal to one in all time
  periods for individuals that are eventually treated

- xformla:

  A optional one sided formula for additional covariates that will be
  adjusted for. E.g ~ age + education. Additional covariates can also be
  passed by name using the x paramater.

- t:

  The 3rd time period in the sample. Treated individuals should be
  treated in this time period and untreated individuals should not be
  treated. The code attempts to enforce this condition, but it is good
  try to handle this outside the panel.qtet method.

- tmin1:

  The 2nd time period in the sample. This should be a pre-treatment
  period for all individuals in the sample.

- tmin2:

  The 1st time period in the sample. This should be a pre-treatment
  period for all individuals in the sample.

- tname:

  The name of the column containing the time periods

- data:

  A data.frame containing all the variables used

- idname:

  The individual (cross-sectional unit) id name

- probs:

  A vector of values between 0 and 1 to compute the QTET at

- iters:

  The number of iterations to compute bootstrap standard errors. This is
  only used if se=TRUE

- alp:

  The significance level used for constructing bootstrap confidence
  intervals

- method:

  The method for including covariates, should either be "QR" for
  quantile regression or "pscore" for propensity score

- se:

  Boolean whether or not to compute standard errors

- retEachIter:

  Boolean whether or not to return list of results from each iteration
  of the bootstrap procedure (default is FALSE). This is potentially
  useful for debugging but can cause errors due to running out of
  memory.

- pl:

  Whether or not to compute standard errors in parallel

- cores:

  Number of cores to use if computing in parallel

- biters:

  Number of bootstrap iterations; alias for `iters` matching the
  `did`/`ptetools` naming convention. If both are supplied, `biters`
  takes precedence.

- cl:

  Number of cores for parallel bootstrap; alias for `pl`/`cores`.
  `cl = 1` (default) runs sequentially; `cl > 1` enables parallelism.

## Value

`QTE` object

## References

Callaway, Brantly and Tong Li. “Quantile Treatment Effects in Difference
in Differences Models with Panel Data.” Working Paper, 2019.

## See also

[`panel_qtt`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
for the modern replacement supporting staggered adoption.

## Examples

``` r
# See ?panel_qtt for the modern replacement.
if (FALSE) { # \dontrun{
data(lalonde)

## Run the panel.qtet method on the experimental data with no covariates
pq1 <- panel.qtet(re ~ treat,
  t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
  data = lalonde.exp.panel, idname = "id", se = FALSE,
  probs = seq(0.05, 0.95, 0.05)
)
summary(pq1)

## Run the panel.qtet method on the observational data with no covariates
pq2 <- panel.qtet(re ~ treat,
  t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
  data = lalonde.psid.panel, idname = "id", se = FALSE,
  probs = seq(0.05, 0.95, 0.05)
)
summary(pq2)

## Run the panel.qtet method on the observational data conditioning on
## age, education, black, hispanic, married, and nodegree.
## The propensity score will be estimated using the default logit method.
pq3 <- panel.qtet(re ~ treat,
  t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
  xformla = ~ age + I(age^2) + education + black + hispanic + married + nodegree,
  data = lalonde.psid.panel, idname = "id", se = FALSE, method = "pscore",
  probs = seq(0.05, 0.95, 0.05)
)
summary(pq3)

pq4 <- panel.qtet(re ~ treat,
  t = 1978, tmin1 = 1975, tmin2 = 1974, tname = "year",
  xformla = ~ age + I(age^2) + education + black + hispanic + married + nodegree,
  data = lalonde.psid.panel, idname = "id", se = FALSE, method = "qr",
  probs = seq(0.05, 0.95, 0.05)
)
summary(pq4)
} # }
```
