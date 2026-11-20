# Rearranged conditional quantiles from a quantile regression

Predicted conditional quantiles at each row of `newdata`, with each row
sorted. Sorting is the rearrangement of Chernozhukov, Fernandez-Val and
Galichon (2010): it makes the conditional quantile function monotone
when the quantile regressions cross, and leaves it unchanged when they
don't.

Used instead of `predict(fit, type = "Fhat"/"Qhat", stepfun = TRUE)`,
which fails when `newdata` has a single row and returns a non-monotone
`Fhat` when quantiles cross (issue \#9).

## Usage

``` r
qr_sorted_quantiles(fit, newdata)
```

## Arguments

- fit:

  an `rqs` object (`rq()` with a vector of `tau`).

- newdata:

  data frame of covariate values; if missing, the data used to fit the
  model.

## Value

n x length(fit\$tau) matrix; row i holds the sorted predicted quantiles
for row i of `newdata`.

## References

Chernozhukov, V., I. Fernandez-Val and A. Galichon (2010). “Quantile and
Probability Curves Without Crossing.” Econometrica 78(3), 1093-1125.
