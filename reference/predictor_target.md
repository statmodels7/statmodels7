# The Response on the Scale of a Predictor

Returns the response carried onto the scale of one equation's linear
predictor.
[`modelterms7::term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.html)
takes it to estimate a term's own parameters from the data, which is the
one thing that a term cannot work out for itself.

## Usage

``` r
predictor_target(spec, p)
```

## Arguments

- spec:

  The specification.

- p:

  The distribution parameter naming the equation.

## Value

A numeric vector, one value per observation, or `NULL`.

## Details

A target exists only where the response reads the parameter directly,
which `params_interpretation` records as a mean or a location. A scale
or a shape has no per-observation reading, so `NULL` is returned.

The response is passed on the scale of the predictor, \\g(y)\\, because
a term that estimates its parameters from the raw response would work on
the wrong scale. The contribution of the other terms of the equation is
not subtracted, and the offset of the equation is.

A response that lies on a bound of the domain of the parameter is moved
half way from the bound to the smallest (or, for an upper bound, the
largest) response inside the domain, because a link is not defined at
its bound. For counts under a log link this gives one half when the
smallest positive count is 1.

## See also

[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md),
[`modelterms7::term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.html)
