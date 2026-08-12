# Which Terms Rewrite the Likelihood

The parameter, name and shape of every structural term in a
specification.

## Usage

``` r
statmod_structural(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A list of entries with `param`, `term` and `kind`, possibly empty.

## Details

The shape is `"filter"` for a term that produces a predictor, which the
layer adds to the equation the term sits in, and `"loglik"` for one
whose contribution is a likelihood mixed over states and cannot be
written as a predictor at all.

## See also

[`statmod_filter_at`](https://statmodels7.github.io/statmodels7/reference/statmod_filter_at.md)
