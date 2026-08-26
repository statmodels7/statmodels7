# The Tolerance a Criterion Asks For

Returns the largest tolerance a criterion asks for, descending into a
combined rule and taking the maximum over its parts.

## Usage

``` r
criterion_tol(crit)
```

## Arguments

- crit:

  An optimizers7 criterion, possibly a combined one built by
  `crit_any()` or `crit_all()`.

## Value

A single number.
[`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html)'s
own default where the criterion carries no tolerance at all, read from
that function's formals.

## Details

The maximum is the right reduction because a combined rule stops when
**any** part fires, so the loop it drives is no more precise than its
loosest term. A caller of
[`method_budget()`](https://statmodels7.github.io/statmodels7/reference/method_budget.md)
wants the precision the loop will actually deliver.

## See also

[`method_budget()`](https://statmodels7.github.io/statmodels7/reference/method_budget.md),
its caller.
