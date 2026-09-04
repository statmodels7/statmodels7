# What the Inner Method Says About How to Fit

The information matrix, its approximation and the budget, read off the
inner method in one place.

## Usage

``` r
inner_settings(method)
```

## Arguments

- method:

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

## Value

A list of four: `expected` (a logical), `approx` (a string), `maxit` and
`tol`.

## Details

Every route that fits the coefficients reads these, and reading them in
one place keeps a refit inside a path or a fold on the same terms as the
fit the caller asked for. Hard-coding them instead ran cross-validation
at a tolerance a hundred times tighter than
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)'s
own, which cost 26 per cent in time and answered a question the caller
had not asked.

An optimizers7 optimizer gets the OBSERVED information. It minimizes
\\-\ell + \rho\\ and
[`optimizers7::minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.html)
documents `he` as that function's second derivative, which is what the
observed information is; the expected one is a different matrix, so a
method asked for a Newton step was performing Fisher scoring under
another name. Where a family writes neither out the two also differ in
cost by orders: the route taken before was the exact sum over the
support, measured at 211 s per call on a 12096-cell Poisson-inverse
gaussian regression against 0.041 s for the observed, so `newton()`
spent an hour building matrices and never finished. `approx` is `"opg"`
for such a method, which is what
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
itself defaults to; it is read only where `expected` is `TRUE`, so for
an optimizer it records a default rather than a choice.

## See also

[`method_budget()`](https://statmodels7.github.io/statmodels7/reference/method_budget.md)
for the last two,
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for where the first two are set.
