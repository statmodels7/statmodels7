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

A single number.

## Details

optimizers7 defaults to 30, which suits an objective costing
microseconds. Every trial here is a penalized refit, and a line search
that is going to fail spends its whole budget finding that out.

The saving is smaller than the evaluation count implies. Every trial
warm-starts from the last accepted point, so as the step shrinks the
refit begins at nearly its own answer: measured, removing 22 evaluations
of 38 removed 2.8 seconds of 30.8, which is 0.13 seconds each against an
average evaluation's 0.81.

What it costs is nothing on the shapes where nothing was wrong. Measured
at 30, 12 and 8 backtracks on a smooth, on two smooths with a random
effect and on a random intercept, all three are unchanged in
evaluations, in criterion to six decimals, in effective degrees of
freedom and in the convergence flag: the resolution
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
supplies stops the search before the budget is ever reached. A
hierarchical break-point model goes from 31 evaluations and 25.6 seconds
to 13 and 20.3, with the criterion 1.3e-04 better and the same degrees
of freedom.

The value is 12 rather than 8 because 12 takes 5.3 seconds of the 6.6
there are to take, and because the criterion 8 gives up is larger: swept
with the optimizer named, so that the resolution does not mask it, the
criterion lost against a budget of 30 is 4.8e-05 at 12, 7.8e-04 at 8 and
3.3e-03 at 6, the last being past the 1.6e-03 that criterion can
resolve.

An optimizer the caller named keeps its own budget, as it keeps its own
stopping rule.

## See also

[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
