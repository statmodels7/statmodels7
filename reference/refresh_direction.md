# What One Direction Costs a Moving Block

The two quantities every consumer of \\\partial X/\partial\beta\\ reads
for a given direction: the block's own derivative along it, and the
third derivative already contracted in it and carried onto the term's
columns.

## Usage

``` r
refresh_direction(spec, design, M, params, npar, offs, d3, tv, v, units)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- M:

  The matrix the trace is taken against.

- params, npar, offs:

  The block bookkeeping.

- d3:

  The third derivatives on the link scale.

- tv:

  The predictors of the direction.

- v:

  The direction, over the stacked coefficients.

- units:

  The refreshable terms, from
  [`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

## Value

A list, one entry per unit, with `D` and `A`.

## Details

Both depend on one direction alone, never on a pair, so they are built
once per hyperparameter and combined in the Hessian's pair loop:
computing them inside that loop would repeat an \\O(np^2)\\ product for
every pair.

## See also

[`contract3_refresh()`](https://statmodels7.github.io/statmodels7/reference/contract3_refresh.md),
[`trace_refresh4()`](https://statmodels7.github.io/statmodels7/reference/trace_refresh4.md)
