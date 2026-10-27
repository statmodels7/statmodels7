# The Settings Read From the Inner Method

Reads the information matrix, its approximation and the budget from the
inner method in one place.

## Usage

``` r
inner_settings(method, distrib = NULL)
```

## Arguments

- method:

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- distrib:

  The distribution, read only to settle an
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  left at
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  through
  [`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md);
  such a method with no distribution is an error rather than a guess.

## Value

A list of four: `expected` (a logical), `emethod` (an object), `maxit`
and `tol`.

## Details

Every route that fits the coefficients reads these settings, and reading
them in one place keeps a refit inside a path or a fold on the same
terms as the fit of the caller.

An optimizers7 optimizer receives the observed information. It minimizes
\\-\ell + \rho\\ and
[`optimizers7::minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.html)
documents `he` as that function's second derivative, which is the
observed information; the expected information is a different matrix, so
a Newton step on it would be Fisher scoring under another name. Where a
family writes neither out, the expected information is also far more
expensive to form than the observed one. `emethod` is `NULL` for such a
method, the default of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
itself; it is read only where `expected` is `TRUE`, so for an optimizer
it records a default and not a choice.

## See also

[`method_budget()`](https://statmodels7.github.io/statmodels7/reference/method_budget.md)
for the last two,
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for where the first two are set.
