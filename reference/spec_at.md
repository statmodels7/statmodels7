# Rebuild a Specification Against New Data

Returns the fit's own specification when `data` is `NULL`, and one
carrying the fitted terms reapplied to the new rows otherwise. This is
the one place the two cases are told apart, so every route that accepts
a `data` argument treats `NULL` alike.

## Usage

``` r
spec_at(fit, data, need_response = TRUE)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- data:

  A data frame, or `NULL` for the data the model was fitted to.

- need_response:

  `TRUE` where the response must be present, as for a log-likelihood;
  `FALSE` for a prediction, new data routinely having no response
  column.

## Value

A
[`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md):
`fit@spec` itself when `data` is `NULL`, and
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)'s
result otherwise.
