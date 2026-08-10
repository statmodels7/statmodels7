# Rebuild a Specification Against New Data

Returns the fit's specification when `data` is `NULL`, and one built
against the new data otherwise.

## Usage

``` r
spec_at(fit, data, need_response = TRUE)
```

## Arguments

- fit:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- data:

  A data frame or `NULL`.

- need_response:

  Whether the response must be there. A likelihood needs it; a
  prediction does not, and new data routinely has no response column.

## Value

A
[`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).
