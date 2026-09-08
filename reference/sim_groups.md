# The Groups a Printed Truth Collapses

One entry per scalar parameter and one per family of names that share a
stem, in the order the names arrive.

## Usage

``` r
sim_groups(nm)
```

## Arguments

- nm:

  The names.

## Value

A list of entries, each with `label` and `idx`.

## Details

A structural term's parameters are `omega`, `alpha1` and the like where
nothing is developed, and `omega.(Intercept)`, `omega.random.1` and so
on where something is. The first are read one by one and the second are
a block, so the printed form is grouped at the first dot, the scalars
coming out as groups of one and needing no special case.

## See also

[`print.StatmodSim()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodSim.md)
