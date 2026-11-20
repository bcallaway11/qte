# Conditional CDF and quantile function from sorted QR predictions

Evaluate the conditional CDF (`qr_cdf`) or quantile function
(`qr_quantile`) for each row of `q_sorted` by linear interpolation
between the estimated quantiles `(tau_k, Q(tau_k))`. Interpolating,
rather than using a step function, avoids rounding ranks down to the
quantile grid, which biases CiC-type counterfactuals.

Beyond the grid the functions are flat: `qr_cdf` returns `tau_1` below
the first quantile and `tau_K` above the last, and `qr_quantile` returns
the first (last) quantile for `p` below `tau_1` (above `tau_K`). In
simulations, extrapolating with exponential tails reduced bias only
slightly under normal errors and badly hurt heavy-tailed (t3) designs,
and estimating extra extreme quantiles added more variance than it
removed bias.

## Usage

``` r
qr_cdf(q_sorted, taus, y)

qr_quantile(q_sorted, taus, p)
```

## Arguments

- q_sorted:

  matrix from
  [`qr_sorted_quantiles()`](https://bcallaway11.github.io/qte/reference/qr_sorted_quantiles.md).

- taus:

  quantile levels, `fit$tau`.

- y, p:

  numeric vectors with one value per row of `q_sorted`.

## Value

numeric vector, one value per row of `q_sorted`.
