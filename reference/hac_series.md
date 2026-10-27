# The Series of a HAC Variance

Returns the rows of each independent series, in the order of the series.

## Usage

``` r
hac_series(fit, design, sw)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

- sw:

  The
  [HacSandwich](https://statmodels7.github.io/statmodels7/reference/HacSandwich.md)
  object.

## Value

A list of integer vectors of row numbers.

## Details

A structural term gives them in its blueprint. For a model without one
they are the rows of each level of `cluster`, ordered by `time`, with
ties kept in the order of the data.
