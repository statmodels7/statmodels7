# The Derivative Row of an Equation Carrying a Likelihood Mixed Over States

The derivative of the posterior-weighted predictor of the equation that
a structural term of the likelihood shape (a regime term or a marginal
break-point term) shifts, one row per observation, in the coefficients
of every equation and the term's free parameters on their unconstrained
scale.

## Usage

``` r
loglik_shift_jacobian(spec, design, coef)
```

## Arguments

- spec:

  The specification.

- design:

  Its design.

- coef:

  The coefficients.

## Value

A list with `X` (the derivative rows) and `key` (the names of their
columns in the variance matrix), restricted to the free columns, or
`NULL` where the model carries no such term.

## Details

The predictor is the static part plus the posterior mean of the term's
shift. Its derivative is the equation's design row plus the derivative
of that mean, which
[`modelterms7::term_shift_jacobian()`](https://statmodels7.github.io/modelterms7/reference/term_shift_jacobian.html)
returns: the posterior moves with every coordinate, through the
emissions of every equation and through the term's own parameters, and
the shifts move with the term's own parameters. The family's score at
each component's predictor is supplied here, as for
[`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md).

## See also

[`structural_se_columns()`](https://statmodels7.github.io/statmodels7/reference/structural_se_columns.md),
[`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md)
