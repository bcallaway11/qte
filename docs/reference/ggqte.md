# ggqte

**Deprecated.** Use `autoplot()` on a `QTE` object instead.

## Usage

``` r
ggqte(
  qteobj,
  main = "",
  ylab = "QTE",
  ylim = NULL,
  ybreaks = NULL,
  xbreaks = c(0.1, 0.3, 0.5, 0.7, 0.9),
  setype = "pointwise",
  alp = qteobj$alp
)
```

## Arguments

- qteobj:

  a QTE object

- main:

  optional title

- ylab:

  optional y axis label

- ylim:

  optional limits of y axis

- ybreaks:

  optional breaks in y axis

- xbreaks:

  optional breaks in x axis

- setype:

  options are "pointwise", "uniform" or both

- alp:

  significance level override

## Value

a ggplot object
