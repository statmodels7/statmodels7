# Every Penalty in a Model, Without the Design

Returns the same enumeration as
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
without the column positions, for the callers that need to know which
penalties exist before a design has been built.

## Usage

``` r
statmod_penalty_keys(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A list of entries with `param`, `term`, `key`, `within` (positions among
the term's own parameters), `penalty`, `fixed`, `n_values`, `values`,
`min_ratio`, `search` and `ids`. The entry of a covariance class also
carries `class` and `params`.

## See also

[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
