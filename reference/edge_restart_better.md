# Whether a Restart From the Edge Improved the Fit

Compares a fit restarted by
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
after
[`edge_violations()`](https://statmodels7.github.io/statmodels7/reference/edge_violations.md)
with the one it replaces: on the criterion that chose the
hyperparameters where one ran, in its own direction, and on the
penalized objective otherwise.

## Usage

``` r
edge_restart_better(new, old, crit, outer_criterion, sparse_criterion)
```

## Arguments

- new, old:

  The two results of
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
  or
  [`statmod_select()`](https://statmodels7.github.io/statmodels7/reference/statmod_select.md).

- crit:

  The criterion `old` reached, `NA` where none ran.

- outer_criterion, sparse_criterion:

  As passed to
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md).

## Value

`TRUE` where `new` is strictly better.
