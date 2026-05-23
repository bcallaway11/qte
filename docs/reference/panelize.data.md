# panelize.data

get data in correct format for using panel methods in qte package

## Usage

``` r
panelize.data(data, idname, tname, t, tmin1, tmin2 = NULL)
```

## Arguments

- data:

  A data.frame containing all the variables used

- idname:

  The individual (cross-sectional unit) id name

- tname:

  The name of the column containing the time periods

- t:

  The 3rd time period in the sample. Treated individuals should be
  treated in this time period and untreated individuals should not be
  treated. The code attempts to enforce this condition, but it is good
  try to handle this outside the panel.qtet method.

- tmin1:

  The 2nd time period in the sample. This should be a pre-treatment
  period for all individuals in the sample.

- tmin2:

  The 1st time period in the sample. This should be a pre-treatment
  period for all individuals in the sample.

## Value

data.frame
