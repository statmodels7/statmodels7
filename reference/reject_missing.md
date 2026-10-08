# Reject Missing Values in the Data a Fit Reads

Signals an error naming the response, or each variable of the data that
an equation names, where it holds a missing value, with the rows. A fit
does not drop rows: in a time series or a panel a dropped row changes
what the neighbouring rows mean, so the choice of which rows to remove
is the caller's.

## Usage

``` r
reject_missing(response, equations, data)
```

## Arguments

- response:

  The response, as evaluated.

- equations:

  The equations, one formula per parameter.

- data:

  The data frame.

## Value

`NULL`, invisibly, when nothing is missing.
