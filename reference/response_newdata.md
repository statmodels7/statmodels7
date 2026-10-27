# New Data for an Effect on the Response Scale

Completes a grid with every variable of the model, those of the terms
and those of the offsets: the variables in `at` take the given values,
and the others their median or first level in `data`.

## Usage

``` r
response_newdata(fit, data, grid, at)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- data:

  The data frame of the fit.

- grid:

  The grid of the effect.

- at:

  A named list, or `NULL`.

## Value

A data frame with the rows of `grid`.
