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
    #>  1794.343        653.17   514.1535    3074.533 *
    #>
    #>
    #> QTE:
    #>  Tau      QTE Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1    0.000     0.0000         0.0000      0.0000
    #>  0.2    0.000   127.4163      -347.6936    347.6936
    #>  0.3  929.884   373.4064       -89.0671   1948.8350
    #>  0.4 1176.944   971.9214     -1475.2345   3829.1217
    #>  0.5 1091.468   870.9200     -1285.0976   3468.0334
    #>  0.6 1466.690   906.4408     -1006.8042   3940.1852
    #>  0.7 1810.696  1022.2255      -978.7522   4600.1437
    #>  0.8 2300.454  1165.8815      -881.0025   5481.9100
    #>  0.9 2856.068  1882.8732     -2281.9142   7994.0509
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
    #>  -12537.37      1966.364  -16391.37   -8683.368 *
    #>
    #>
    #> QTE:
    #>  Tau        QTE Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1      0.000   703.0556      -2398.343   2398.3433
    #>  0.2  -6404.679  1833.9759     -12660.946   -148.4118 *
    #>  0.3  -8523.773  2256.9373     -16222.894   -824.6517 *
    #>  0.4 -11114.085  2678.9947     -20252.977  -1975.1939 *
    #>  0.5 -12684.516  2911.2438     -22615.681  -2753.3497 *
    #>  0.6 -13532.837  2715.4178     -22795.979  -4269.6946 *
    #>  0.7 -14181.889  2066.6976     -21232.043  -7131.7352 *
    #>  0.8 -18285.313  2180.2109     -25722.696 -10847.9295 *
    #>  0.9 -22580.785  3887.0913     -35840.872  -9320.6968 *
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
    #>  -4543.922        853.17  -6216.105    -2871.74 *
    #>
    #>
    #> QTT:
    #>  Tau         QTT Std. Error [ 95% Simult.  Conf. Band]
    #>  0.1      0.0001   128.0586      -343.4861    343.4862
    #>  0.2  -1002.7420   820.4547     -3203.4130   1197.9290
    #>  0.3  -3400.5673  1572.4851     -7618.3778    817.2432
    #>  0.4  -5009.2491  1223.5366     -8291.0895  -1727.4086 *
    #>  0.5  -4602.4652  1115.6798     -7595.0062  -1609.9243 *
    #>  0.6  -5229.1454  1269.7220     -8634.8672  -1823.4237 *
    #>  0.7  -5507.4720  1415.9424     -9305.3946  -1709.5493 *
    #>  0.8  -6885.7529  1466.8178    -10820.1365  -2951.3693 *
    #>  0.9 -10517.0625  2503.9534    -17233.3109  -3800.8141 *
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
