# Whether a Restart From the Edge Improved the Fit

Compares a fit restarted by
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
(after
[`edge_violations()`](https://statmodels7.github.io/statmodels7/reference/edge_violations.md),
or from a score-driven term's own start) with the one it replaces. The
comparison is made on the criterion that chose the hyperparameters,
where one ran, in the direction of that criterion, and on the penalized
objective otherwise.

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
  The element `criterion_by` of `old`, the
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  whose value is `criterion`, gives the criterion.

- crit:

  The criterion `old` reached, `NA` where none ran.

- outer_criterion, sparse_criterion:

  As passed to
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
  read only where `old` does not record the criterion that produced
  `crit`; the outer criterion is then used if there is one, and the
  sparse criterion otherwise.

## Value

`TRUE` where `new` is strictly better.
