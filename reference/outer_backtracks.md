# The Outer Line Search's Backtracking Budget

How many backtracks a line search inside
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
may take before it gives up, where the optimizer is the one this package
chose.

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
applies only to the line search of an optimizer that this package chose;
the resolution that
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
supplies ends the search before the budget on most shapes.

An optimizer that the caller named keeps its own budget, as it keeps its
own stopping rule.

## See also

[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
