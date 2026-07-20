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

The Quantile Treatment Effect at level $\tau$ is
$QTE(\tau) = Q_{Y(1)}(\tau) - Q_{Y(0)}(\tau)$, where $Q_{Y(d)}(\tau)$ is
the $\tau$-th quantile of the potential outcome distribution under
treatment status $d$. Mapping $QTE(\tau)$ across $\tau \in (0,1)$ traces
how the treatment shifts the entire outcome distribution. The Quantile
Treatment Effect on the Treated ($QTT$) conditions on the treated group.

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
  biters = 50
)
summary(res_exp)
#> 
#> Overall ATE:  
#>       ATE    Std. Error     [ 95%  Conf. Int.]  
#>  1794.343      655.1211   510.3294    3078.357 *
#> 
#> 
#> QTE:
#>  Tau      QTE Std. Error [ 95% Simult.  Conf. Band]  
#>  0.1    0.000     0.0000         0.0000      0.0000  
#>  0.2    0.000   127.5450      -249.9836    249.9836  
#>  0.3  929.884   361.5139       221.3297   1638.4382 *
#>  0.4 1176.944   950.4180      -685.8414   3039.7286  
#>  0.5 1091.468   961.2442      -792.5361   2975.4718  
#>  0.6 1466.690   956.3520      -407.7249   3341.1059  
#>  0.7 1810.696  1047.2267      -241.8308   3863.2223  
#>  0.8 2300.454  1161.5595        23.8390   4577.0686 *
#>  0.9 2856.068  1811.4848      -694.3767   6406.5134  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

``` r
autoplot(res_exp, ylab = "QTE (earnings, 1978)")
```

![QTE curve under random
assignment](intro_files/figure-html/random-qte-plot-1.png)

The confidence band indicates that the QTE is positive across much of
the distribution, with the effect somewhat larger in the lower
quantiles.

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
- `"or"` — outcome regression (conditional quantile regression)
- `"aipw"` — augmented IPW (doubly robust; recommended)

``` r
xf <- ~ age + I(age^2) + education + black + hispanic + married + nodegree

res_psid <- unc_qte(
  yname      = "re78",
  dname      = "treat",
  data       = lalonde.psid,
  xformla    = xf,
  est_method = "aipw",
  target     = "qte",
  probs      = seq(0.1, 0.9, 0.1),
  biters     = 100
)
summary(res_psid)
#> 
#> Overall ATE:  
#>        ATE    Std. Error     [ 95%  Conf. Int.]  
#>  -13194.78      1560.465  -16253.24   -10136.33 *
#> 
#> 
#> QTE:
#>  Tau        QTE Std. Error [ 95% Simult.  Conf. Band]  
#>  0.1      0.000   588.7718      -1153.972    1153.972  
#>  0.2  -6404.679  1583.9763      -9509.216   -3300.142 *
#>  0.3  -8523.773  1953.7527     -12353.058   -4694.488 *
#>  0.4 -11114.085  2294.7224     -15611.659   -6616.512 *
#>  0.5 -12684.516  2724.8872     -18025.196   -7343.835 *
#>  0.6 -13532.837  2583.7614     -18596.916   -8468.757 *
#>  0.7 -14181.889  1779.0101     -17668.685  -10695.093 *
#>  0.8 -18285.313  1578.1361     -21378.403  -15192.223 *
#>  0.9 -22580.785  3257.9287     -28966.208  -16195.362 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

``` r
autoplot(res_psid, ylab = "QTE (earnings, 1978)")
```

![QTE curve under
unconfoundedness](intro_files/figure-html/obs-qte-plot-1.png)

## QTT under unconfoundedness

Setting `target = "qtt"` estimates the Quantile Treatment Effect on the
Treated — the distributional effect for the subpopulation that actually
received treatment.

``` r
res_qtt <- unc_qte(
  yname      = "re78",
  dname      = "treat",
  data       = lalonde.psid,
  xformla    = xf,
  est_method = "aipw",
  target     = "qtt",
  probs      = seq(0.1, 0.9, 0.1),
  biters     = 100
)
summary(res_qtt)
#> 
#> Overall ATT:  
#>        ATT    Std. Error     [ 95%  Conf. Int.]  
#>  -4685.583       891.021  -6431.952   -2939.214 *
#> 
#> 
#> QTT:
#>  Tau         QTT Std. Error [ 95% Simult.  Conf. Band]  
#>  0.1      0.0001   109.2933      -214.2109    214.2110  
#>  0.2  -1002.7420   849.1294     -2667.0051    661.5211  
#>  0.3  -3400.5673  1676.6309     -6686.7035   -114.4311 *
#>  0.4  -5009.2491  1294.4558     -7546.3359  -2472.1623 *
#>  0.5  -4602.4652  1059.1306     -6678.3231  -2526.6073 *
#>  0.6  -5229.1454  1358.6095     -7891.9711  -2566.3198 *
#>  0.7  -5507.4720  1633.4368     -8708.9492  -2305.9947 *
#>  0.8  -6885.7529  1652.9272    -10125.4307  -3646.0752 *
#>  0.9 -10517.0625  2308.9602    -15042.5414  -5991.5836 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
```

``` r
autoplot(res_qtt, ylab = "QTT (earnings, 1978)")
```

![QTT curve under
unconfoundedness](intro_files/figure-html/qtt-plot-1.png)

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

Firpo, Sergio. 2007. “Efficient Semiparametric Estimation of Quantile
Treatment Effects.” *Econometrica* 75 (1): 259–76.
