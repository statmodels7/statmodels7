# The Outer Line Search's Backtracking Budget

Returns the number of backtracks that a line search inside
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
may take before it stops, where the optimizer is the one that
[`outer_default_optimizer()`](https://statmodels7.github.io/statmodels7/reference/outer_default_optimizer.md)
selects.

## Usage

``` r
outer_backtracks()
```

## Value

A single number, 12.

## Details

optimizers7 defaults to 30, which suits an objective that costs
microseconds. Every trial here is a penalized refit, and a line search
that fails spends its whole budget before it stops.

Every trial warm-starts from the incumbent, the best usable point found
so far, so as the step shrinks the refit begins near its own answer and
a rejected backtrack costs less than an average evaluation. The budget
applies only to the line search of an optimizer selected by
[`outer_default_optimizer()`](https://statmodels7.github.io/statmodels7/reference/outer_default_optimizer.md);
the resolution that
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
supplies ends the search before the budget on most shapes.

An optimizer that the caller named keeps its own budget, as it keeps its
own stopping rule.

## See also

[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
