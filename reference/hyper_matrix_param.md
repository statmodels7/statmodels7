# The Matrix Parameter of a Multivariate Prior

Returns the chart, the names of its free values among the family's
parameters, and whether the matrix is a precision, read from a
multivariate family under any number of
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
wrappers.

## Usage

``` r
hyper_matrix_param(d)
```

## Arguments

- d:

  A distributions7 object.

## Value

A list with `chart`, `names` and `inverted`, or `NULL` when `d` is not a
multivariate family that carries one matrix parameter.
