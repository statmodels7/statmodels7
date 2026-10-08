# A Least-Squares Start for a Location on the Identity Link

Fits the parametric columns of the first parameter's equation to the
response by least squares (weighted by the prior weights, net of that
equation's offset), and fits the distribution without covariates to the
residuals shifted by the fitted intercept. Returns the coefficients of
those columns and the intercepts of every parameter, or `NULL` where the
start does not apply.

## Usage

``` r
location_ols_start(spec, design)
```

## Arguments

- spec:

  The specification.

- design:

  The design.

## Value

`NULL`, or a list with `param` (the first parameter's name), `cols`
(positions in its coefficient vector), `coef` (their starting values)
and `eta0` (link-scale intercepts, as from
[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)).

## Details

It applies when the family is univariate, its first parameter is a
location or a mean on the identity link, that parameter's equation has
an intercept and at least one other parametric column, and every term of
the equation is a
[modelterms7::LinparTerm](https://statmodels7.github.io/modelterms7/reference/LinparTerm.html),
a
[modelterms7::SmoothTerm](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.html),
a
[modelterms7::RandomTerm](https://statmodels7.github.io/modelterms7/reference/RandomTerm.html)
or a
[modelterms7::PenalizedTerm](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.html).
The penalized blocks start at zero, as they do from the intercept-only
start, and a term with a start of its own (a break-point, a nonlinear
term, a filter) leaves the equation on the intercept-only start. An
aliased column takes the coefficient zero.

## See also

[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md),
[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
