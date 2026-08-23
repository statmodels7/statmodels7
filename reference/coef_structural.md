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

  The fit.

- p:

  The distribution parameter.

- readable:

  Whether to report the quantities or the coordinates.

## Value

A named numeric vector, empty where the equation carries no structural
term.

## Details

A structural term rewrites the likelihood rather than adding columns to
a design, so its parameters are in no block and were in no reading of
[`coef`](https://rdrr.io/r/stats/coef.html): a model whose whole
predictor is a score-driven filter answered with an empty vector. They
are named from the term's label as every other coefficient of the term
is.

A held parameter is one an intercept in the same equation carries and is
not estimated; it is reported under either reading, at the value it is
held at, because leaving it out would make the vector shorter than the
term's own parameter count.
