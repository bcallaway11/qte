README
================
Brantly Callaway
2023-09-16

<!-- README.md is generated from README.Rmd. Please edit that file -->

# qte <img src="man/figures/logo.png" align="right" height="139" alt="" />

The R `qte` package implements many methods used, especially in
economics, to estimate quantile treatment effects. These include the
case where treatment is randomly assigned, under selection on
observables, under a Difference in Differences Assumtpion.

The package is available on CRAN and can be loaded as follows

``` r
  library(qte)
```

The following example shows how to use the `ci.qte` method in the `qte`
package using data about an experimental job training program.

``` r
data(lalonde)
 jt.cia <- ci.qte(re78 ~ treat,
   xformla=~age + education + black + hispanic + married + nodegree,
   data=lalonde.psid,
   probs=seq(0.05,0.95,0.05), se=T)
 summary(jt.cia)
```

More examples and details about other functions in the package can be
found at the package’s [website](http://bcallaway11.github.io/qte/)
