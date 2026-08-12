# Which Terms Recompute Their Own Block

The parameter and name of every term whose design block is a function of
its own coefficients, in the order the design holds them.

## Usage

``` r
statmod_refreshable(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A list of entries with `param` and `term`, possibly empty.

## See also

[`statmod_design_at`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md),
[`refreshes_own_block`](https://statmodels7.github.io/statmodels7/reference/refreshes_own_block.md)
