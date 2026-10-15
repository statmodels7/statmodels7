# Starting Values From the Intercept-Only Fit

The answer of the default strategy: the intercept-only fit, the
least-squares start of a location where
[`location_ols_start()`](https://statmodels7.github.io/statmodels7/reference/location_ols_start.md)
applies, and the start that each term declares through
[`modelterms7::term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.html).
See
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md).

## Arguments

- strategy:

  A `StartIntercepts` object.

- spec, design, obj, ...:

  As in
  [`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md).

## Value

A named list of numeric vectors.
