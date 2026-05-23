# bootstrap

`bootstrap` is a method that bootstraps standard errors for methods in
the qte package

## Usage

``` r
bootstrap(qteparams, qteest, func)
```

## Arguments

- qteparams:

  QTEparams object containing the parameters passed to ciqte

- qteest:

  an initial estimate of the qte

- func:

  the internal function that computes estimates of the qte, should be
  the same as what was used to produce qteest

## Value

SE object
