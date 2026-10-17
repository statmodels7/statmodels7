# The Groups That a Printed Truth Collapses

Returns one entry per scalar parameter and one per family of names that
share a stem, in the order in which the names arrive.

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

A structural term's parameters are `omega`, `kappa1` and the like where
a parameter is not developed, and `omega.(Intercept)`, `omega.random.1`
and so on where it is. The first are read one by one and the second are
a block, so the printed form is grouped at the first dot. The scalars
then come out as groups of one and need no special case.

## See also

[`print.StatmodSim()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodSim.md)
