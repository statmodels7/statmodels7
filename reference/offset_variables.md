# Variables of the Offsets

Returns the columns of `data` that the offsets written in the formula of
the fit read.

## Usage

``` r
offset_variables(spec, data)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- data:

  The data frame of the fit.

## Value

A character vector, possibly empty.

## Details

The offsets are re-evaluated on new data when the fit predicts, so a
graph on the response scale needs a value for each of their variables.
