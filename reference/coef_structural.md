# What a Structural Term Contributes to the Coefficients

Returns the own parameters of a structural term of one equation, as
quantities or as the coordinates on which they were estimated.

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

  `TRUE` for the quantities that the term declares, `FALSE` for the
  coordinates on the unconstrained scale.

## Value

A named numeric vector of that term's own parameters. The result is
`numeric(0)` where the equation carries no structural term.

## Details

A structural term rewrites the likelihood and adds no columns to a
design, so its parameters sit in no block.
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
appends them to the vector of the equation under both readings, named
from the term's label as every other coefficient of that term is.

A held parameter is one that an intercept in the same equation carries,
so it is not estimated. It is reported under either reading, at the
value at which it is held. Leaving it out would make the vector shorter
than the term's own parameter count, and the readable vector would no
longer correspond to the rows of `vcov(readable = TRUE)`, where the held
level has a row of missing values.

## See also

[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md),
the caller,
[`statmod_latent()`](https://statmodels7.github.io/statmodels7/reference/statmod_latent.md)
for a latent-state term's smoothed states.
