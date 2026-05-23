# compute.unc_qte

Computes the QTE or QTT under unconfoundedness. Dispatches to IPW, OR
(outcome regression), or AIPW based on `qp$est_method` and `qp$target`.
Called internally by `unc_qte` and during each bootstrap iteration.

## Usage

``` r
compute.unc_qte(qp)
```

## Arguments

- qp:

  QTEparams object

## Value

QTE object
