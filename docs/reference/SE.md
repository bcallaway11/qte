# SE

Class for Standard Error Objects

## Usage

``` r
SE(
  qte.se = NULL,
  ate.se = NULL,
  qte.upper = NULL,
  qte.lower = NULL,
  ate.upper = NULL,
  ate.lower = NULL,
  c = NULL,
  probs = NULL
)
```

## Arguments

- qte.se:

  The QTE Standard Error

- ate.se:

  The ATE Standard Error

- qte.upper:

  The QTE upper CI

- qte.lower:

  The QTE lower CI

- ate.upper:

  The ATE upper CI

- ate.lower:

  The ATE lower CI

- c:

  The critical value from a KS-type statistic used for creating uniform
  confidence bands

- probs:

  The values at which the QTE is computed
