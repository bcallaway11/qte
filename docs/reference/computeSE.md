# computeSE

Computes standard errors from bootstrap results. This function is called
by several functions in the qte package

## Usage

``` r
computeSE(bootIters, qteobj, alp = 0.05)
```

## Arguments

- bootIters:

  List of bootstrap iterations

- alp:

  The significance level used for constructing bootstrap confidence
  intervals

## Value

SEObj
