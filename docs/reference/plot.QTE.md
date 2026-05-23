# plot.QTE

Plots a `QTE` object using
[`autoplot.QTE`](https://bcallaway11.github.io/qte/reference/autoplot.QTE.md).

## Usage

``` r
# S3 method for class 'QTE'
plot(x, cband = TRUE, ylab = "QTE", ...)
```

## Arguments

- x:

  a `QTE` object, as returned by
  [`unc_qte`](https://bcallaway11.github.io/qte/reference/unc_qte.md).

- cband:

  logical; if `TRUE` (default), show the uniform confidence band. If
  `FALSE`, show pointwise intervals.

- ylab:

  label for the y-axis. Default `"QTE"`.

- ...:

  passed to
  [`autoplot.QTE`](https://bcallaway11.github.io/qte/reference/autoplot.QTE.md).

## Value

invisibly returns the `ggplot` object.
