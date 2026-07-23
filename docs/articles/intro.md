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
    #>  1794.343      614.3984   590.1444    2998.542 *
    #>
    #>
    #> QTE:
    #>  Tau      QTE Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1    0.000     0.0000         0.0000      0.0000
    #>  0.2    0.000   158.4803      -310.6157    310.6157
    #>  0.3  929.884   400.0082       145.8823   1713.8856 *
    #>  0.4 1176.944   943.3150      -671.9198   3025.8070
    #>  0.5 1091.468   831.0748      -537.4088   2720.3446
    #>  0.6 1466.690   842.6293      -184.8326   3118.2136
    #>  0.7 1810.696  1032.9801      -213.9081   3835.2995
    #>  0.8 2300.454  1054.5595       233.5552   4367.3524 *
    #>  0.9 2856.068  1795.8551      -663.7429   6375.8796
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

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
```

    #>
    #> Overall ATE:
    #>        ATE    Std. Error     [ 95%  Conf. Int.]
    #>  -13194.78       1669.31  -16466.57   -9922.995 *
    #>
    #>
    #> QTE:
    #>  Tau        QTE Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1      0.000   636.5657      -1247.646    1247.646
    #>  0.2  -6404.679  1637.5230      -9614.165   -3195.193 *
    #>  0.3  -8523.773  2123.4289     -12685.617   -4361.929 *
    #>  0.4 -11114.085  2398.9413     -15815.924   -6412.247 *
    #>  0.5 -12684.516  2756.7181     -18087.584   -7281.448 *
    #>  0.6 -13532.837  2739.1129     -18901.399   -8164.274 *
    #>  0.7 -14181.889  2214.9031     -18523.019   -9840.759 *
    #>  0.8 -18285.313  2239.9340     -22675.503  -13895.123 *
    #>  0.9 -22580.785  3682.3739     -29798.105  -15363.464 *
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

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
```

    #>
    #> Overall ATT:
    #>        ATT    Std. Error     [ 95%  Conf. Int.]
    #>  -4685.583      835.1563  -6322.459   -3048.707 *
    #>
    #>
    #> QTT:
    #>  Tau         QTT Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1      0.0001     9.4440       -18.5099     18.5100
    #>  0.2  -1002.7420   705.1096     -2384.7314    379.2474
    #>  0.3  -3400.5673  1505.4521     -6351.1992   -449.9354 *
    #>  0.4  -5009.2491  1072.6921     -7111.6870  -2906.8111 *
    #>  0.5  -4602.4652   936.3507     -6437.6789  -2767.2516 *
    #>  0.6  -5229.1454  1282.4542     -7742.7095  -2715.5814 *
    #>  0.7  -5507.4720  1268.5601     -7993.8040  -3021.1399 *
    #>  0.8  -6885.7529  1501.3663     -9828.3768  -3943.1290 *
    #>  0.9 -10517.0625  2368.3837    -15159.0092  -5875.1159 *
    #> ---
    #> Signif. codes: `*' confidence band does not cover 0

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
