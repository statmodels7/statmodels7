# Where an Equation's Intercept Is

The position of the intercept of one equation's parametric block, or
`NA` where it has none.

## Usage

``` r
parametric_intercept(spec, design, p)
```

## Arguments

- spec:

  The specification.

- design:

  The design.

- p:

  The distribution parameter naming the equation.

## Value

An integer position into the equation's coefficient vector, or
`NA_integer_`.

## Details

The position searched for is a column of the parametric block, which is
a stricter test than a coefficient whose name ends in `(Intercept)`.
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
names the intercept of each of its own parameters the same way, and
those coefficients lie on the charts of the parameters of the term and
not on the scale of the predictor, so a model written `y ~ 0 + nl(...)`
has one of them first. Writing the intercept-only fit there would set
the parameter to the inverse link of a value that was never on its
scale, and for a parameter on a log link the starting value would be far
outside the range of the data.

## See also

[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md),
[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
