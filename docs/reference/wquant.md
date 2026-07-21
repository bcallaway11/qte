# Weighted type-1 quantile

Smallest `y` with weighted CDF `>= p`. Uses a small numerical tolerance
on the comparison so that a `p` landing exactly on a floating-point
rounding boundary (common with uniform weights, e.g. `1/n`) picks the
intended order statistic instead of flipping based on last-bit
differences in summation across platforms.

## Usage

``` r
wquant(y_vec, w_vec, probs)
```

## Arguments

- y_vec:

  numeric vector of values.

- w_vec:

  numeric vector of weights, same length as `y_vec`.

- probs:

  numeric vector of probabilities in `[0, 1]`.

## Value

numeric vector, same length as `probs`.
