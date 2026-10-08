# The Coordinates of the Nested Random Effects Set Aside

For each row of
[`nested_random_terms()`](https://statmodels7.github.io/statmodels7/reference/nested_random_terms.md),
the positions of the sub-term's coefficients in its equation's design,
and those of the first level.

## Usage

``` r
nested_columns(spec, design, rows)
```

## Arguments

- spec, design:

  The specification and its design.

- rows:

  Rows of
  [`nested_random_terms()`](https://statmodels7.github.io/statmodels7/reference/nested_random_terms.md).

## Value

A list, one entry per row, with `param`, `cols` and `first`.
