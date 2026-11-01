# The Driving Score's Derivatives as the Fifth Order Reads Them

Returns the arrays that
[`modelterms7::term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.html)
reads for a filter in one equation: the derivatives of the driving score
of orders one to five in the predictors of every equation, one row per
observation, the filter's own parameter included.

## Usage

``` r
jfp_fifth_arrays(params, ap, Vs, dr, n)
```

## Arguments

- params:

  The distribution parameters.

- ap:

  The index of the filter's parameter.

- Vs:

  The static rows, one matrix per parameter.

- dr:

  What
  [`filter_driving()`](https://statmodels7.github.io/statmodels7/reference/filter_driving.md)
  returns, to the sixth order.

- n:

  The number of observations.

## Value

A list with `H`, `D3`, `D4`, `D5`, `D6` (matrices with one column per
tuple of parameters, the first index running fastest), `Vs` and `ap`.
