# The Budget and the Stopping Rule of the Alternation

Reads the iteration budget and the tolerance off the method that fits
the smooth block, which is where they are set.

## Usage

``` r
method_budget(method)
```

## Arguments

- method:

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

## Value

A list of two: `maxit`, the budget in passes, and `tol`, the relative
tolerance on the objective.

## Details

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
carries neither a `maxit` nor a `tol` of its own. An argument accepted
and ignored is worse than one that signals an error, and a second copy
would be exactly that: a caller setting `iwls(maxit = 20)` alongside a
loose `maxit = 100` would be obeyed in one and never told about the
other. distributions7's `fit_distrib()` shed the same pair for the same
reason.

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
carries both directly. An optimizers7 optimizer carries `maxit` and a
`criterion`, and the tolerance is taken as the **largest** one in the
criterion tree: a combined rule stops at whichever of its parts fires
first, so the alternation should not ask for more precision than the
loop inside it can deliver.

A criterion carrying no tolerance at all leaves the default of
[`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html),
read from that function's own formals and never copied as a number, so a
change there reaches here.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`criterion_tol()`](https://statmodels7.github.io/statmodels7/reference/criterion_tol.md)
for the tree walk.
