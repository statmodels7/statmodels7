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
has no `maxit` and no `tol`. Both are read from `inner_optimizer`, so
that one setting governs the inner fit and the alternation.

[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
carries both directly. An optimizers7 optimizer carries `maxit` and a
`criterion`, and the tolerance is read from the criterion tree by
[`criterion_tol()`](https://statmodels7.github.io/statmodels7/reference/criterion_tol.md):
the largest tolerance of the parts of a rule built by `crit_any()`,
which stops when any part holds, and the smallest of the parts of a rule
built by `crit_all()`, which stops only when every part holds. The
alternation then requires no more precision than the loop inside it
delivers.

A criterion carrying no tolerance at all takes the default of
[`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html),
read from that function's own formals and not copied as a number, so
that a change there reaches this function.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`criterion_tol()`](https://statmodels7.github.io/statmodels7/reference/criterion_tol.md)
for the tree walk.
