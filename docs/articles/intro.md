# Quantile Treatment Effects in R

## What are Quantile Treatment Effects?

The Average Treatment Effect (ATE) or Average Treatment Effect on the
Treated (ATT) summarizes the causal effect of a treatment by a single
number: the mean difference in outcomes between treated and untreated
units. This is often a natural quantity, but it can miss important
heterogeneity. If a job training program substantially raises earnings
at the bottom of the distribution while having little effect at the top,
the ATE obscures this. If a minimum wage increase compresses the wage
distribution, an average effect conceals the compression.

The Quantile Treatment Effect at level $`\tau`$ is
$`QTE(\tau) = Q_{Y(1)}(\tau) -
Q_{Y(0)}(\tau)`$, where $`Q_{Y(d)}(\tau)`$ is the $`\tau`$-th quantile
of the potential outcome distribution under treatment status $`d`$.
Mapping $`QTE(\tau)`$ across $`\tau \in (0,1)`$ traces how the treatment
shifts the entire outcome distribution. The Quantile Treatment Effect on
the Treated ($`QTT`$) conditions on the treated group.

This vignette demonstrates
[`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md),
the `qte` package’s estimator for cross-sectional settings (no panel
data required). For DiD-based estimators and staggered treatment
adoption, see
[`vignette("panel-estimators")`](https://bcallaway11.github.io/qte/articles/panel-estimators.md)
and
[`vignette("staggered-adoption")`](https://bcallaway11.github.io/qte/articles/staggered-adoption.md).

``` r
library(qte)
library(ggplot2)
set.seed(42)
data(lalonde)
xf <- ~ age + I(age^2) + education + black + hispanic + married + nodegree
```

## Random assignment

Under random assignment, no covariate adjustment is needed: the
unconditional quantiles of the treated and control outcome distributions
identify the QTE directly.
[`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
implements the semiparametrically efficient estimator of Firpo
([2007](#ref-firpo-2007)). We use the experimental Lalonde data
(`lalonde.exp`), where treatment was randomly assigned.

``` r
res_exp <- unc_qte(
  yname  = "re78",
  dname  = "treat",
  data   = lalonde.exp,
  target = "qte",
  probs  = seq(0.1, 0.9, 0.1),
  biters = 100
)
summary(res_exp)
```

    #>
    #> Overall ATE:
    #>       ATE    Std. Error     [ 95%  Conf. Int.]
    #>  1794.343      787.3904   251.0863      3337.6 *
    #>
    #>
    #> QTE:
    #>  Tau      QTE Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1    0.000     0.0000         0.0000      0.0000
    #>  0.2    0.000   160.9787      -496.5387    496.5387
    #>  0.3  929.884   400.1089      -304.2519   2164.0199
    #>  0.4 1176.944   961.4600     -1788.6798   4142.5670
    #>  0.5 1091.468   914.5023     -1729.3145   3912.2503
    #>  0.6 1466.690   921.5603     -1375.8621   4309.2431
    #>  0.7 1810.696  1091.5274     -1556.1205   5177.5120
    #>  0.8 2300.454  1380.8276     -1958.7089   6559.6165
    #>  0.9 2856.068  2189.5705     -3897.6618   9609.7985
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0
    #> Note: the outcome looks discrete (308 distinct values; 31% at 0); quantile effects may be unreliable.

Earnings are zero for 31% of this sample, so
[`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
warns that the outcome looks discrete, and
[`summary()`](https://rdrr.io/r/base/summary.html) repeats the note.
Quantiles in the mass point are zero in both groups, so the QTE there is
zero by construction; interpret the low quantiles with care.

``` r
autoplot(res_exp, ylab = "QTE (earnings, 1978)")
```

![QTE curve under random
assignment](intro_files/figure-html/random-qte-plot-1.png)

The point estimates are zero at the lowest quantiles, inside the mass
point, and then increase with $`\tau`$, from about \$900 at the 30th
percentile to about \$2,900 at the 90th. With 445 participants, though,
the uniform band covers zero at every quantile, even though the ATE of
about \$1,800 is statistically significant: a band that holds
simultaneously over nine quantiles is much wider than a single interval.

## Unconfoundedness

When treatment is not randomly assigned, the unconfoundedness assumption
(conditional independence of potential outcomes given covariates)
provides identification ([Firpo 2007](#ref-firpo-2007)). We use the
observational Lalonde data (`lalonde.psid`), which combines the treated
group from the experiment with a comparison group from the PSID.

[`unc_qte()`](https://bcallaway11.github.io/qte/reference/unc_qte.md)
supports three estimation methods:

- `"ipw"` — inverse probability weighting (propensity score
  re-weighting)
- `"reg"` — outcome regression (conditional quantile regression)
- `"aipw"` — augmented IPW (doubly robust; recommended)

``` r
xf <- ~ age + I(age^2) + education + black + hispanic + married + nodegree

res_psid <- unc_qte(
  yname      = "re78",
  dname      = "treat",
  data       = lalonde.psid,
  xformula   = xf,
  est_method = "aipw",
  target     = "qte",
  probs      = seq(0.1, 0.9, 0.1),
  biters     = 100
)
summary(res_psid)
```

    #>
    #> Overall ATE:
    #>        ATE    Std. Error     [ 95%  Conf. Int.]
    #>  -12535.52      1713.692   -15894.3   -9176.749 *
    #>
    #>
    #> QTE:
    #>  Tau        QTE Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1      0.000   661.1683      -2461.476   2461.4758
    #>  0.2  -6402.223  1714.8165     -12786.345    -18.1013 *
    #>  0.3  -8523.451  1937.3801     -15736.160  -1310.7431 *
    #>  0.4 -11114.190  2334.3133     -19804.648  -2423.7316 *
    #>  0.5 -12682.669  2639.0571     -22507.663  -2857.6741 *
    #>  0.6 -13533.232  2569.7020     -23100.024  -3966.4411 *
    #>  0.7 -14181.863  1930.2224     -21367.924  -6995.8024 *
    #>  0.8 -18285.823  1844.5771     -25153.033 -11418.6125 *
    #>  0.9 -22580.784  4489.0837     -39293.277  -5868.2909 *
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0
    #> Note: the outcome looks discrete (740 distinct values; 12% at 0); quantile effects may be unreliable.

``` r
autoplot(res_psid, ylab = "QTE (earnings, 1978)")
```

![QTE curve under
unconfoundedness](intro_files/figure-html/obs-qte-plot-1.png)

These estimates are large and negative, far from the experimental
benchmark above. The QTE is about the whole population, which here is
mostly PSID households with much higher earnings than the program’s
participants, and these covariates (which leave out earnings before the
program) do not make the two groups comparable ([LaLonde
1986](#ref-lalonde-1986); [Dehejia and Wahba
1999](#ref-dehejia-wahba-1999)).

## QTT under unconfoundedness

Setting `target = "qtt"` estimates the Quantile Treatment Effect on the
Treated — the distributional effect for the subpopulation that actually
received treatment.

``` r
res_qtt <- unc_qte(
  yname      = "re78",
  dname      = "treat",
  data       = lalonde.psid,
  xformula   = xf,
  est_method = "aipw",
  target     = "qtt",
  probs      = seq(0.1, 0.9, 0.1),
  biters     = 100
)
summary(res_qtt)
```

    #>
    #> Overall ATT:
    #>        ATT    Std. Error     [ 95%  Conf. Int.]
    #>  -4543.928        836.89  -6184.202   -2903.654 *
    #>
    #>
    #> QTT:
    #>  Tau         QTT Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1      0.0001    91.9581      -233.1124    233.1126
    #>  0.2  -1002.7420   759.8658     -2928.9911    923.5071
    #>  0.3  -3400.6253  1603.5495     -7465.6003    664.3498
    #>  0.4  -5009.2481  1010.3860     -7570.5624  -2447.9339 *
    #>  0.5  -4602.5048   805.2498     -6643.8014  -2561.2081 *
    #>  0.6  -5229.1200  1109.1159     -8040.7129  -2417.5271 *
    #>  0.7  -5505.9404  1319.7726     -8851.5450  -2160.3357 *
    #>  0.8  -6885.7327  1417.9261    -10480.1553  -3291.3102 *
    #>  0.9 -10517.0895  2256.7599    -16237.9437  -4796.2353 *
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0
    #> Note: the outcome looks discrete (740 distinct values; 12% at 0); quantile effects may be unreliable.

``` r
autoplot(res_qtt, ylab = "QTT (earnings, 1978)")
```

![QTT curve under
unconfoundedness](intro_files/figure-html/qtt-plot-1.png)

The QTT only needs a counterfactual for the participants, so it is much
smaller in magnitude than the QTE, but it is still negative and grows
with $`\tau`$, and the band excludes zero from the 40th percentile up.
It remains well below the experimental estimates, a reminder that
unconfoundedness given these covariates is a strong assumption for this
comparison group.

## Reading the output

[`summary()`](https://rdrr.io/r/base/summary.html) prints the overall
ATT, the QTE/QTT at each quantile, and standard errors with confidence
intervals.
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
returns a `ggplot` object so you can add layers:

``` r
autoplot(res_qtt) +
  ggplot2::labs(title = "QTT — Lalonde (observational)",
                subtitle = "AIPW with pre-treatment covariates")
```

By default
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
shows a uniform confidence band (`cband = TRUE`), which provides
simultaneous coverage over all quantile levels — a stronger guarantee
than pointwise intervals. Pass `cband = FALSE` for pointwise intervals
instead.

Standard errors are computed via the empirical bootstrap. The `biters`
argument controls the number of iterations (default 100). Parallel
computation is available via the `cl` argument.

## References

Dehejia, Rajeev, and Sadek Wahba. 1999. “Causal Effects in
Nonexperimental Studies: Reevaluating the Evaluation of Training
Programs.” *Journal of the American Statistical Association* 94 (448):
1053–62.

Firpo, Sergio. 2007. “Efficient Semiparametric Estimation of Quantile
Treatment Effects.” *Econometrica* 75 (1): 259–76.

LaLonde, Robert. 1986. “Evaluating the Econometric Evaluations of
Training Programs with Experimental Data.” *American Economic Review*,
604–20.
