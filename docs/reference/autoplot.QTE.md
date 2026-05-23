# autoplot.QTE

Plot a `QTE` object as a quantile treatment effect curve with optional
confidence bands.

## Usage

``` r
# S3 method for class 'QTE'
autoplot(object, cband = TRUE, ylab = "QTE", ...)
```

## Arguments

- object:

  a `QTE` object, as returned by
  [`unc_qte`](https://bcallaway11.github.io/qte/reference/unc_qte.md).

- cband:

  logical; if `TRUE` (default), show the uniform confidence band stored
  in `object$qte.upper` / `object$qte.lower`. If `FALSE`, show pointwise
  intervals computed from `object$qte.se`.

- ylab:

  label for the y-axis. Default `"QTE"`.

- ...:

  unused.

## Value

a `ggplot` object.
