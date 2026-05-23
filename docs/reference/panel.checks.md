# panel.checks

Validates that the panel data structure is compatible with the panel
estimators in the qte package. Specifically: checks that the treated
group composition is constant across all time periods, warns if
covariates appear to vary over time, and stops with an informative
message if the treated group is empty or changes.

## Usage

``` r
panel.checks(qp)
```

## Arguments

- qp:

  QTEparams object

## Value

Invisibly returns `NULL`. Called for its side effects (warnings and
errors).
