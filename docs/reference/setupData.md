# setupData

`setupData` extracts treated/untreated subsets, ECDFs, and auxiliary
vectors from a `QTEparams` object and **injects them into the caller's
environment** via `assign(..., envir=parent.frame())`. This
side-effectful pattern is intentional — all `compute.*` functions rely
on the variables created here being available in their local scope.

## Usage

``` r
setupData(qteParams)
```

## Arguments

- qteParams:

  A `QTEparams` object holding the function parameters.

## Value

Invisibly returns `NULL`. The function's effect is the assignment of
data objects into the calling environment.
