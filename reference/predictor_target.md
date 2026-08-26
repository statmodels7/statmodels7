# The Response on the Scale of a Predictor

The response carried onto the scale of one equation's linear predictor.
[`modelterms7::term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.html)
takes it to estimate a term's own parameters from the data, which is the
one thing a term cannot work out for itself.

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

It exists only where the response reads the parameter directly, which
`params_interpretation` says: a mean or a location. A scale or a shape
has no per-observation reading, so `NULL` is returned and nothing is
invented for the occasion.

The scale matters and is not a detail. Measured on a Poisson whose
predictor is a logistic growth curve with \\\phi = 4\\, a term handed
the raw response estimated \\\phi\\ between 52.7 and 54.8 over five
samples, and one handed \\g(y)\\ estimated it between 4.03 and 4.07.
What does not matter, measured on the same shape beside another term, is
the other terms' contribution: subtracting it moved the estimate from
39.68 to 39.87 against a truth of 40, so the response is passed as it
stands and nothing is residualized.

A response sitting on a bound of the parameter's own domain is moved
half way to the nearest admissible value, a link not being defined at
its bound. For counts under a log link that is the classical half.

## See also

[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md),
[`modelterms7::term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.html)
