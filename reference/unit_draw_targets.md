# Where a Penalty's Coordinates Sit in the Vectors a Simulation Draws

One target per coordinate the penalty covers, in the penalty's own
order: which vector it belongs to, which parameter's coefficients, and
the position within them.

## Usage

``` r
unit_draw_targets(u, design, params)
```

## Arguments

- u:

  One entry of
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md).

- design:

  The design.

- params:

  The distribution parameter names.

## Value

A list with `space`, `param` and `pos`, or `NULL` where the entry cannot
be addressed.

## Details

A penalty over a structural term's own parameters is read among those
parameters, on the unconstrained scale, and a penalty over an equation's
coefficients among those coefficients. A covariance class spans several
members and its penalty reads them interleaved group by group, which is
what its `index` already records, so the stacked positions are turned
back into a parameter and a column with the offsets the design gives.

A class whose halves are in both vectors at once is refused rather than
half-drawn: those coordinates fall back to the plain draw, which is
[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)'s
rule for everything no prior reaches.

## See also

[`rstatmod_truth()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_truth.md)
