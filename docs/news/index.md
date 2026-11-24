# Changelog

## qte 2.1.0 (development version)

### Breaking changes

- Removed the deprecated two-period functions, as announced in 2.0.0:
  `CiC()`, `QDiD()`, `MDiD()`, `ddid2()`, `panel.qtet()`, `ci.qte()`,
  `ci.qtet()` and `ggqte()`. Use
  [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md),
  [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md),
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md),
  [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  (with `target = "qtt"` for the QTT;
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  for staggered adoption) and
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  instead.
- Removed the exported helper `panelize.data()`, which only served those
  functions.
- Removed the exported internal `QTEparams()` parameter container.
  [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  now has self-contained estimation and bootstrap code and no longer
  uses it.
- `formula.tools`, `data.table` and `parallel` are no longer
  dependencies.
- Requires `ptetools (>= 1.1.0)` and `BMisc (>= 1.4.11)`.
- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)’s
  covariate argument `xformla` is renamed `xformula`, matching the other
  estimators. `xformla` still works for now, with a deprecation warning.
  The new `clustervars` argument comes right after `weightsname`, so
  calls passing `probs` or later arguments by position need updating.
- `est_method = "or"` (outcome regression) in
  [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  and
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  is renamed `"reg"`, the name `did` uses, and the name in the new
  [`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md).
  `"or"` still works, with a deprecation warning.

### New features

- [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md),
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  and
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  gain a `clustervars` argument, as in
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html),
  to cluster standard errors at a level above the unit (ptetools
  [\#6](https://github.com/bcallaway11/qte/issues/6)). With the
  empirical bootstrap, whole clusters are resampled.

- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  also gains `clustervars`: its bootstrap then resamples whole clusters
  instead of rows.

- QTT results (`gt_type = "qtt"`) from the panel estimators now also
  carry the ATT, aggregated with the same weights as the QTT, and the
  curve QTT(tau) - ATT with standard errors and a uniform band, from the
  same bootstrap draws (ptetools 1.1.0, ptetools
  [\#27](https://github.com/bcallaway11/qte/issues/27)). Its p-value
  tests QTT(tau) = ATT at every tau: whether the effect is more than a
  location shift. See `autoplot(res, type = "qtt_minus_att")`. This
  includes
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)’s
  default long copula, which has its own aggregation.

- New
  [`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
  (and its workhorse
  [`cdf_did_gt()`](https://bcallaway11.github.io/qte/reference/cdf_did_gt.md)):
  the QTT and ATT under parallel trends in the *distribution functions*
  of untreated potential outcomes (Roth and Sant’Anna 2023), with
  staggered adoption, for panel data or repeated cross sections. The
  counterfactual distribution is computed exactly on the data’s support
  and rearranged into a valid distribution function; the ATT equals the
  usual DiD ATT. Unlike the rank- and quantile-based estimators, it is
  identified for discrete outcomes. Covariates enter through
  `est_method = c("aipw", "ipw", "reg")`: a doubly robust, an
  inverse-probability-weighted, and an outcome-regression (distribution
  regression) version of the conditional assumption (Sant’Anna and Zhao
  2020).

- The estimators now warn, once per call, when the outcome looks
  discrete: more than 10% of observations at one value, or at most 20
  distinct values
  ([\#17](https://github.com/bcallaway11/qte/issues/17)). Quantile
  treatment effects usually require a continuous outcome: quantiles jump
  at mass points, rank-based identification only gives bounds, and
  bootstrap inference can be unreliable. The warning has class
  `qte_discrete_outcome`, the thresholds are the options
  `qte.mass_point_tol` and `qte.discrete_n_unique`, and the diagnostics
  are returned as `outcome_discreteness`.
  [`cic_gt()`](https://bcallaway11.github.io/qte/reference/cic_gt.md) no
  longer warns about “extreme ranks” once per group-time cell and
  bootstrap draw.

- The covariate-adjusted estimators gain `cond_dist = c("qr", "dr")`:
  [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md),
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  (and their `*_gt()` workhorses), and the outcome-regression and AIPW
  branches of
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  and
  [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md).
  With covariates, `"dr"` estimates conditional outcome distributions by
  distribution regression (a logit of 1{Y \<= y} on the covariates at a
  grid of thresholds) instead of quantile regression. It handles
  outcomes with mass points (e.g. many zeros) naturally and avoids
  truncating each unit’s conditional distribution at its 1st and 99th
  percentiles. The default stays `"qr"`, and the `"qr"` results are
  unchanged.

### Bug fixes (these change results)

- Quantile regressions (`cond_dist = "qr"`, the default) now use the
  Frisch-Newton interior point method (`quantreg::rq(method = "fn")`)
  instead of the Barrodale-Roberts simplex. With a mass point in the
  outcome (e.g. earnings with many zeros), the simplex could cycle
  forever and hang R, uninterruptibly. Estimates are unchanged up to
  solver tolerance when the quantile regression solution is unique.

- [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) and
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  could pick the wrong order statistic for a weighted quantile when a
  CDF value landed exactly on a probability (e.g. 35 weights of 1/500
  summing to 0.07) and summation rounding put it a hair below. Which one
  was picked depended on the platform (results differed between macOS on
  Apple Silicon and Linux) and, in
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md),
  on rounding even within one platform. They now use the same small
  tolerance as the package’s other weighted quantiles.
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  estimates can change slightly.

- When a treated unit’s covariates aren’t spanned by the comparison
  group (e.g. a dummy for a region with treated units only), the
  conditional outcome models extrapolate. With `cond_dist = "dr"` this
  happened silently: the dummy got a zero coefficient, so those treated
  units got the reference category’s conditional distribution. With
  `cond_dist = "qr"` the estimators stopped with an uninformative
  `Singular design matrix`. Both now drop such columns and warn once per
  call, naming the covariates, the group-time cells and the (weighted)
  share of treated units affected
  ([\#16](https://github.com/bcallaway11/qte/issues/16)). This applies
  to every estimator with covariates, including
  [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md),
  [`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
  and
  [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md).
  Dropped columns that don’t matter (a factor level no unit in the cell
  has, a duplicated covariate) give no warning, and covariates with full
  support give the same results as before. Use
  [`ptetools::check_overlap()`](https://rdrr.io/pkg/ptetools/man/check_overlap.html)
  for overlap diagnostics by cell.

- Factor covariates (and
  [`factor()`](https://rdrr.io/r/base/factor.html) or character terms in
  `xformula`) with a level that the comparison group lacks in some cell
  no longer stop the covariate-adjusted estimators with a “new levels”
  error.

- With panel data and covariates that change over time, the estimators
  now condition on each unit’s covariates from the base period of each
  group-time cell, as `ptetools` and `did` do. Before,
  [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) and
  [`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md)
  fit the comparison group’s post-period distribution on post-period
  covariates, and
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  conditioned on covariates measured after treatment.
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  now fixes covariates at the earliest period each cell uses; it
  previously mixed two pre-treatment periods. Results with
  time-invariant covariates, and with repeated cross sections, are
  unchanged.

- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)’s
  uniform confidence band could be absurdly wide (by orders of
  magnitude) when a requested quantile fell inside a mass point of the
  outcome, such as zero earnings. The critical value standardized each
  quantile’s bootstrap draws by an IQR-based scale, which is nearly 0
  there, while the band used the standard deviation. The draws are now
  standardized by the same standard deviation the band uses, and
  quantiles whose draws never vary are left out of the critical value,
  as the panel estimators do (ptetools 1.1.0, ptetools
  [\#33](https://github.com/bcallaway11/qte/issues/33)). Standard errors
  are unchanged; uniform bands change.

- With more than one treated cohort, the aggregated QTTs (overall, event
  study and group) mixed the wrong group-time cells, including
  pre-treatment ones. The aggregation weights were in a different order
  from the cells they were applied to. Fixed in ptetools 1.1.0, which
  qte now requires. This affected
  [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md),
  [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md),
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  and `panel_qtt(pre_copula = "short")` with `gt_type = "qtt"` or
  `"qott"`; ATTs were not affected, and neither was
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)’s
  default `pre_copula = "long"`, whose aggregation is in qte. On the
  package’s test data the median QTT moves by up to 0.13, toward the
  true effect.

- Fixed covariate-adjusted
  [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) and
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  (and their `*_gt()` workhorses)
  ([\#9](https://github.com/bcallaway11/qte/issues/9)):

  - They no longer fail when a group has a single treated unit.
  - Conditional quantiles are now rearranged (sorted) before use, so the
    conditional CDFs are monotone when quantile regressions cross.
  - Conditional ranks are no longer one quantile-grid step (0.01) too
    low. This came from quantreg’s
    `predict(..., type = "Fhat", stepfun = TRUE)` and biased the
    counterfactual distribution downward.
  - Conditional CDFs and quantile functions are now linearly
    interpolated between the estimated quantiles, not step functions, so
    ranks are no longer rounded down to the quantile grid. Beyond the
    first and last estimated quantiles they stay flat (no tail
    extrapolation).

  Covariate-adjusted estimates change as a result. In a simulation where
  conditional CiC holds exactly, the bias of the ATT falls from 0.15 to
  0.02. Estimates without covariates are unchanged.

- [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) /
  [`ddid_gt()`](https://bcallaway11.github.io/qte/reference/ddid_gt.md)
  with covariates now integrate the counterfactual over the *treated*
  group’s covariate distribution, by reweighting control units with the
  propensity-score odds p(X)/(1 - p(X))
  ([\#10](https://github.com/bcallaway11/qte/issues/10)). Previously the
  controls’ covariate distribution was used, which biased the ATT and
  QTT whenever treatment depends on the covariates (bias 0.84 in a
  simple simulated design; now about 0). Results without covariates are
  unchanged.

- `est_method = "aipw"` in
  [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  and
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  now reports a doubly-robust ATT/ATE: the mean of the AIPW
  counterfactual distribution, using the outcome model’s conditional
  mean ([\#11](https://github.com/bcallaway11/qte/issues/11)).
  Previously the reported ATT/ATE was the IPW estimate, so it ignored
  the outcome model and wasn’t doubly robust. The AIPW quantile curves
  are unchanged.

- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)‘s
  `cband = TRUE` (the default) now reports a uniform confidence band,
  `qte` +/- `c` \* `qte.se`, built the same way as the panel estimators’
  uniform bands. Previously `cband` was ignored: `qte.lower` /
  `qte.upper` were always pointwise intervals, though
  [`summary()`](https://rdrr.io/r/base/summary.html) labelled them
  “Simult.”. The bootstrap critical value `c` was already computed and
  stored. Bands are typically wider under the default; use
  `cband = FALSE` for pointwise intervals.
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  label and draw the band correctly, including for `QTE` objects saved
  by earlier versions.

- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  with `weightsname`: each bootstrap draw now resamples units together
  with their weights. Previously the weights stayed in their original
  row order while the rows were resampled, so standard errors and
  confidence bands with weights were wrong. Point estimates are
  unchanged.

- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  now uses `weightsname` throughout: in the propensity-score model, in
  the outcome-regression fits, and in the AIPW distribution functions
  (previously only in the final averaging).
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)’s
  propensity-score model is now weighted too. With integer weights, the
  estimates now equal those from the unweighted fit on row-replicated
  data. Results without weights are unchanged.

- With groups and time periods coded as calendar years, naming the group
  column `G` or the time column `period` gave `NaN` estimates (ptetools
  [\#13](https://github.com/bcallaway11/qte/issues/13)). Fixed in
  ptetools, which now keeps its internal columns separate from the
  user’s;
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)’s
  subset function uses those internal columns. Requires
  `ptetools (>= 1.1.0)`. Results for other column names are unchanged.

- [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  now treats `anticipation` in units of the time variable, as
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)
  and ptetools do: a cell is post-treatment from `g - anticipation` on,
  pre1 is the last period at least `anticipation` before `g`, and
  not-yet-treated comparison units must not be treated or anticipating
  treatment by the cell’s latest period. It previously counted
  anticipation in periods and ignored it for the comparison group.
  Results without anticipation are unchanged.

- `panel_qtt(gt_type = "qtt")` with the default `pre_copula = "long"` no
  longer fails (“no valid post-treatment cells”) or drops valid cells
  when time periods are unequally spaced, e.g. `lalonde.psid.panel`
  (1974, 1975, 1978). The aggregation step checked whether the earlier
  base period exists in calendar years, while the estimator chooses base
  periods by position in the panel. Results with equally spaced periods
  are unchanged.

- With `weightsname`, the panel estimators’ overall and event-study
  results (ATT and QTT) now combine cohorts by their weighted shares, as
  `did` does. They used unweighted counts of units, so results were
  wrong whenever the average weight differed across cohorts (ptetools
  [\#23](https://github.com/bcallaway11/qte/issues/23); for
  `panel_qtt(pre_copula = "long")`, which aggregates in qte,
  [\#13](https://github.com/bcallaway11/qte/issues/13)). Integer weights
  now give the same estimates as duplicating units. Results without
  weights are unchanged.

- [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) with
  covariates no longer fails with “object ‘pre_ctrl’ not found”. The
  outcome regressions looked up their weights where
  [`lm()`](https://rdrr.io/r/stats/lm.html) could not see them, so
  [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md)
  errored for any non-empty `xformula`.

### Minor improvements

- `msm` is no longer a suggested package; nothing used it.
- Aggregated QTTs are now computed exactly: the group-time cells’ CDFs
  are mixed from every support point and the mixture inverted exactly,
  instead of on a 1,000-point grid of the pooled outcome
  ([\#19](https://github.com/bcallaway11/qte/issues/19); ptetools
  [\#25](https://github.com/bcallaway11/qte/issues/25) for the
  estimators that use ptetools’ aggregation, and
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)’s
  long-copula aggregation in qte). QTT estimates change slightly (by a
  few thousandths on the package’s tests and `mpdta`), and weighted data
  now gives the same QTTs as the equivalent row-replicated data.
  Requires `BMisc (>= 1.4.11)`.
- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  now checks up front that `dname` is a column of `data` coded 0/1, with
  an informative error. A treatment coded e.g. 1/2 previously failed
  inside the propensity-score fit with “y values must be 0 \<= y \<= 1”
  ([\#1](https://github.com/bcallaway11/qte/issues/1)).
- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  drops rows with missing values in the outcome, treatment, weights or
  covariates (with a warning) instead of failing. Its internals were
  rewritten as standalone code; results without weights are unchanged.
- [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) and
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md) with
  covariates and `weightsname` (panel data): the ATT paired the
  counterfactual outcomes, sorted by unit id, with weights in data row
  order ([\#14](https://github.com/bcallaway11/qte/issues/14)). ptetools
  sorts panel data by id, so results were not affected in practice; the
  ATT now uses the same id-ordered weights as the counterfactual
  distribution.
- Invalid sampling weights in the panel estimators (missing, negative,
  or varying over time within a unit) now give a clear error, from
  ptetools ([\#24](https://github.com/bcallaway11/qte/issues/24)),
  instead of `NaN` or silently inconsistent estimates.
- Repeated cross sections
  ([\#18](https://github.com/bcallaway11/qte/issues/18)): the help pages
  and README now say which estimators accept `panel = FALSE`
  ([`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md),
  [`cdf_did()`](https://bcallaway11.github.io/qte/reference/cdf_did.md),
  and
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  with `lagged_outcome_cov = FALSE`) and why
  [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) and
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  need panel data. New tests check that, on panel data, `panel = FALSE`
  gives the panel estimates, and that the estimators work on genuine
  repeated cross sections.
- New tests check that integer weights match row replication for every
  panel estimator, and that rescaling the weights changes nothing
  ([\#15](https://github.com/bcallaway11/qte/issues/15)).
- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)’s
  `pscore.reg` no longer carries copies of the data. The fit kept the
  estimator’s working data twice over, through `$data`/`$model` and
  through its formula’s environment, so a saved `QTE` object for
  `lalonde.psid` shrinks from about 3.6 MB to 0.4 MB.
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  `predict(newdata = )` are unchanged;
  [`anova()`](https://rdrr.io/r/stats/anova.html) and
  [`update()`](https://rdrr.io/r/stats/update.html), which refit the
  model, no longer work on it.

### Performance

- [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md),
  [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md),
  [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md),
  [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md),
  [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  and
  [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  are faster, with the same results:
  - Weighted ranks and quantiles in the `*_gt()` functions use one sort
    instead of a pass over the data for each unit, so they are no longer
    quadratic in the group size.
  - The covariate versions of
    [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) and
    [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md) fit
    the post-period quantile regression for the comparison group once
    per period, not once per cohort, because it is the same for every
    cohort already treated by that period.
  - Only the columns the estimators use are passed to
    [`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html).

  Combined with the faster empirical bootstrap in ptetools 1.1.0,
  unconditional CiC on a 3,000-unit, 30-period panel runs more than 20
  times faster.

## qte 2.0.0

CRAN release: 2026-07-23

This is a major update that extends all estimators to support staggered
treatment adoption. The API has changed substantially — please read this
section before upgrading.

### New functions

- [`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
  — cross-sectional QTE and QTT under unconfoundedness. Supports IPW,
  outcome regression, and doubly robust (AIPW) estimation. Uses a
  `yname`/`dname` interface consistent with the `did`/`ptetools`
  convention. Replaces `ci.qte()` and `ci.qtet()`, which are retained as
  deprecated wrappers.

- [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) /
  [`cic_gt()`](https://bcallaway11.github.io/qte/reference/cic_gt.md) —
  Change in Changes estimator for staggered treatment adoption (Athey
  and Imbens 2006). Supports panel and repeated cross sections data, and
  returns either ATT (`gt_type = "att"`) or a full QTT curve
  (`gt_type = "qtt"`). Replaces `CiC()`, which is retained as a
  deprecated wrapper.

- [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md) /
  [`qdid_gt()`](https://bcallaway11.github.io/qte/reference/qdid_gt.md)
  — Quantile Difference in Differences estimator for staggered treatment
  adoption (Callaway, Li, and Oka 2018). Replaces `QDiD()`, which is
  retained as a deprecated wrapper.

- [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) /
  [`mdid_gt()`](https://bcallaway11.github.io/qte/reference/mdid_gt.md)
  — Mean Difference in Differences estimator for staggered treatment
  adoption (Callaway and Li 2019). Replaces `MDiD()`, which is retained
  as a deprecated wrapper.

- [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) /
  [`ddid_gt()`](https://bcallaway11.github.io/qte/reference/ddid_gt.md)
  — Distributional Difference in Differences estimator for staggered
  treatment adoption (Callaway and Li 2019). Replaces `ddid2()`, which
  is retained as a deprecated wrapper.

- [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md)
  /
  [`panel_qtt_gt()`](https://bcallaway11.github.io/qte/reference/panel_qtt_gt.md)
  — Panel QTT estimator using copula transfer across pre-treatment
  periods (Callaway and Li 2019), with staggered treatment adoption
  support. The `pre_copula` argument controls whether a rolling
  (`"long"`, default) or fixed (`"short"`) base period is used to
  transfer the copula. Replaces `panel.qtet()`, which is retained as a
  deprecated wrapper.

- [`lou_qtt()`](https://bcallaway11.github.io/qte/reference/lou_qtt.md)
  — Lagged-outcome unconfoundedness QTT estimator for staggered
  treatment adoption. Uses the lagged outcome as a control variable to
  identify the QTT. New; no prior equivalent in the package.

### Deprecated functions

The following functions emit a runtime deprecation warning via
[`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) and will be
removed in v2.1 (target: 2026-08). Each points to its modern replacement
in the warning message.

| Deprecated | Replacement |
|----|----|
| `ci.qte()` | `unc_qte(target = "qte")` |
| `ci.qtet()` | `unc_qte(target = "qtt")` |
| `CiC()` | [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) |
| `QDiD()` | [`qdid()`](https://bcallaway11.github.io/qte/reference/qdid.md) |
| `MDiD()` | [`mdid()`](https://bcallaway11.github.io/qte/reference/mdid.md) |
| `ddid2()` | [`ddid()`](https://bcallaway11.github.io/qte/reference/ddid.md) |
| `panel.qtet()` | [`panel_qtt()`](https://bcallaway11.github.io/qte/reference/panel_qtt.md) |

### Removed functions

- `computeDiffSE()` — was dead code; never called anywhere in the
  package. Design notes for a future QTT-vs-ATT test are in
  `dev/NOTES.md`.

- `DR()`, `spatt()`, `bounds()` — removed in earlier 1.5.x development
  work.

- `qtes2mat()`, `diffQ()` — formatting helpers tied to the legacy `QTE`
  class with no internal callers.

### Infrastructure

- Full `testthat` test suite added (83 fast tests; 10 slow bootstrap
  tests skipped by default, enabled via `R_SLOW_TESTS=true`). The
  package previously had no tests.

- All new estimators are backed by
  [`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html). The
  `attgt_fun` hook pattern makes it straightforward to add new
  identification strategies without duplicating bootstrap or aggregation
  code.

- Estimators now pass an explicit `aggte_fun` to
  [`ptetools::pte()`](https://rdrr.io/pkg/ptetools/man/pte.html),
  replacing the earlier implicit dispatch based on `gt_type`.

- `ptetools` added to `Imports` (was missing from `DESCRIPTION` despite
  being used).

- `Hmisc` and `texreg` removed from `Imports` (were unused after earlier
  cleanup).

- License updated from GPL-2 to GPL-3.

------------------------------------------------------------------------

## qte 1.4.0

- Added `cic2()` for Change in Changes with multiple periods and
  variation in treatment timing, using `ptetools` as a backend. This
  function was experimental; it has been superseded by
  [`cic()`](https://bcallaway11.github.io/qte/reference/cic.md) in
  v2.0.0.

## qte 1.3.3

- Bug fix: critical values for uniform confidence bands now correctly
  reflect the `alp` argument rather than being stuck at `alp = 0.05`.

## qte 1.3.2

- Significance level is now passed through in the results.
- Added ability to change significance level in `ggqte()` for pointwise
  confidence intervals without re-estimating.

## qte 1.3.1

CRAN release: 2022-09-01

- Updated maintainer contact information.
- Improvements and updates to vignettes.

## qte 1.3.0

CRAN release: 2019-06-10

- Added first-step quantile regression support for `panel.qtet()`,
  `QDiD()`, and `CiC()`.
- Added checks that input data is in the correct format, particularly
  for panel data.

## qte 1.2.3

- Now using empirical quantiles in the CiC method.
- Added code for computing uniform confidence bands.
- General improvements to the codebase, including combining shared code
  across methods.

## qte 1.2.2

CRAN release: 2018-04-30

- Bug fixes to the CiC method.

## qte 1.2

- Added `ddid2()` to compute quantile treatment effects under a
  Difference in Differences assumption following Callaway and Li (2016).
- Covariates for all methods now specified using formulas (`xformla`
  argument) rather than passing variable names directly.
- Improvements to bootstrapping standard errors.
- Some extensions to computing the ATT via `spatt()`.
- Preliminary version of multi-period ATT moved to the `did` package.
