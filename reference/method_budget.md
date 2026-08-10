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

A list with `maxit` and `tol`.

## Details

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
carries neither a `maxit` nor a `tol` of its own. An argument accepted
and ignored is worse than one that signals an error, and that is what a
second copy would be: a caller setting `iwls(maxit = 20)` and a loose
`maxit = 100` would get one of them with nothing said about the other.
distributions7's `fit_distrib()` shed the same pair for the same reason.

[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
carries both directly. An optimizers7 optimizer carries `maxit` and a
`criterion`; the tolerance is the largest one the criterion tree
contains, since a combined rule stops at whichever of its parts fires
first and the alternation should not ask for more precision than the
loop inside it can deliver. A criterion carrying no tolerance at all
leaves the default of
[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html),
read from that function rather than copied as a number.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
