# Standard Errors of a Linear Map of Coefficients with Missing Variances

Returns \\\sqrt{x_i' V x_i}\\ for each row \\x_i\\ of `X`. A coefficient
with no variance (an entry of `V` that is `NA`) is treated as absent in
the rows that do not read it, and the standard error is `NA` in the rows
that do.

## Usage

``` r
moving_se(X, V)
```

## Arguments

- X:

  A matrix, one row per evaluation.

- V:

  The variance matrix of the coefficients that the columns of `X`
  multiply.

## Value

A numeric vector with one value per row of `X`.
