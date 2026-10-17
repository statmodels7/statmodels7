# Where a Penalty's Coordinates Sit in the Vectors That a Simulation Draws

Returns one target per coordinate that the penalty covers, in the
penalty's own order. Each target holds the vector to which the
coordinate belongs, the parameter whose coefficients contain it and its
position within them.

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
members and its penalty reads them interleaved group by group, in the
order that its `index` already records, so the stacked positions are
turned back into a parameter and a column with the offsets that the
design gives.

A class whose halves are in both vectors at once is not drawn from its
prior: those coordinates fall back to the plain draw, which is
[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)'s
rule for every coordinate that a prior does not reach.

## See also

[`rstatmod_truth()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_truth.md)
