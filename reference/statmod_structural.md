# The Structural Terms of a Specification

Returns the parameter, name and shape of every structural term in a
specification.

## Usage

``` r
statmod_structural(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A list of entries with `param`, `term` and `kind`, possibly empty.

## Details

The shape is `"filter"` for a term that produces a predictor, which the
layer adds to the equation in which the term sits, and `"loglik"` for a
term whose contribution is a likelihood mixed over components (the
states of a regime chain, the intervals or nodes of a marginal
break-point term), which cannot be written as a predictor.

## See also

[`statmod_filter_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_filter_at.md)
