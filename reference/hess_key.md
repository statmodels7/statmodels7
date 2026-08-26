# The Name of a Second-Derivative Component

Locates the \\(a, b)\\ entry of a distribution's Hessian list, which is
keyed by
[`distributions7::hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.html)
and never by position.

## Usage

``` r
hess_key(params, a, b)
```

## Arguments

- params:

  The parameter names, in the family's order, as `distrib@params` gives
  them.

- a, b:

  Indices into `params`, in either order: the Hessian is symmetric and
  this returns the key the list actually holds.

## Value

A single string, one of the names of the list
[`distributions7::distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.html)
returns.

## Details

The key is built by pasting the two parameter names in the order
[`distributions7::hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.html)
uses, and is never parsed back out of a name. That is the discipline
distributions7 records for a family whose parameter name contains the
separator: splitting `"log_scale_mu"` on the underscore gives the wrong
pair, while generating the name from the same enumeration that generated
the list cannot.

## See also

[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md),
its caller.
