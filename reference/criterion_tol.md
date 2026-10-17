# The Tolerance of a Criterion

Returns the tolerance of a criterion, descending into a combined rule:
the largest tolerance of its parts for a rule built by `crit_any()` and
the smallest for a rule built by `crit_all()`.

## Usage

``` r
criterion_tol(crit)
```

## Arguments

- crit:

  An optimizers7 criterion, possibly a combined one built by
  `crit_any()` or `crit_all()`.

## Value

A single number. The value is
[`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html)'s
own default where the criterion carries no tolerance at all, read from
that function's formals.

## Details

A rule built by `crit_any()` stops when any part holds, so the loop that
it drives is no more precise than its loosest part, and the maximum is
taken. A rule built by `crit_all()` stops only when every part holds, so
the loop reaches the precision of its tightest part, and the minimum is
taken. The value is the precision that the loop delivers, and a caller
of
[`method_budget()`](https://statmodels7.github.io/statmodels7/reference/method_budget.md)
needs that value.

## See also

[`method_budget()`](https://statmodels7.github.io/statmodels7/reference/method_budget.md),
its caller.
