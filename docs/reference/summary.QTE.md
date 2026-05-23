# Summary

`summary.QTE` summarizes a `QTE` object, returning formatted data frames
for the overall ATE and the QTE table suitable for printing.

## Usage

``` r
# S3 method for class 'QTE'
summary(object, ...)
```

## Arguments

- object:

  A `QTE` object, as returned by
  [`unc_qte`](https://bcallaway11.github.io/qte/reference/unc_qte.md).

- ...:

  unused.

## Value

A `summary.QTE` object (a list with `overall_ate` and `qte_table` data
frames).
