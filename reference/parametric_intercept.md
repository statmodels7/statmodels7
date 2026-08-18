# Where an Equation's Intercept Is

The position of the intercept of one equation's PARAMETRIC block, or
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

It is a column of the parametric block and not merely a coefficient
whose name ends in `(Intercept)`.
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.html)
names the intercept of each of its own parameters the same way, and
those live on those parameters' charts rather than on the predictor's,
so a model written `y ~ 0 + nl(...)` puts one of them first. Writing the
intercept-only fit there sets the parameter to `linkinv` of a value that
was never on its scale: measured on a logistic growth curve whose
asymptote rides a log link, `mean(y) = 23.9` became a starting \\\phi\\
of `2.5e10`, an objective of `7.0e20` and a gradient of `1.4e21`, on
data whose every scale is ordinary. The lasso path built at those
coefficients then spanned `2.8e15` to `2.8e19` where the block empties
at about 300, so all of its points were the same empty fit and every
subject deviation was estimated as exactly zero.

## See also

[`start_at`](https://statmodels7.github.io/statmodels7/reference/start_at.md),
[`statmod_intercepts`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
