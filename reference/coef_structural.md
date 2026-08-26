# What a Structural Term Contributes to the Coefficients

The own parameters of a structural term of one equation, as quantities
or as the coordinates they were estimated on.

## Usage

``` r
coef_structural(spec, fit, p, readable)
```

## Arguments

- spec:

  The fitted specification.

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- p:

  The distribution parameter naming the equation, a string.

- readable:

  `TRUE` for the quantities the term declares, `FALSE` for the
  coordinates on the unconstrained scale.

## Value

A named numeric vector of that term's own parameters. `numeric(0)` where
the equation carries no structural term.

## Details

A structural term rewrites the likelihood instead of adding columns to a
design, so its parameters sit in no block and appeared in no reading of
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md):
a model whose whole predictor is a score-driven filter answered with an
empty vector. They are named from the term's label, as every other
coefficient of that term is.

A **held** parameter is one an intercept in the same equation carries,
so it is not estimated. It is reported under either reading, at the
value it is held at. Leaving it out would make the vector shorter than
the term's own parameter count and break the correspondence with
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## See also

[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md),
the caller,
[`statmod_latent()`](https://statmodels7.github.io/statmodels7/reference/statmod_latent.md)
for a latent-state term's smoothed states.
