# Does a Structural Term Carry a Penalty of Its Own?

`TRUE` where some penalty of the model covers the parameters of a
structural term rather than a block of design columns.

## Usage

``` r
structural_penalized(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A single logical.

## See also

[`statmod_marginal_full`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
