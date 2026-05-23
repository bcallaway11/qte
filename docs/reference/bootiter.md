# bootiter

`bootiter` that computes a single bootstrap iteration

## Usage

``` r
bootiter(i, qteparams, func)
```

## Arguments

- i:

  the bootstrap iteration, this is not used except for calling by
  apply-type methods

- qteparams:

  QTEparams object containing the parameters passed to ciqte

- func:

  the internal function that computes estimates of the qte, should be
  the same as what was used to produce qteest

## Value

SE object
