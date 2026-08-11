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

A list with `expected`, `approx`, `maxit` and `tol`.

## Details

Every route that fits the coefficients reads these, and reading them in
one place is what keeps a refit inside a path or a fold on the same
terms as the fit the caller asked for. Hard-coding them instead ran
cross-validation at a tolerance a hundred times tighter than
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)'s
own, which cost 26 per cent in time and answered a question the caller
had not asked.
