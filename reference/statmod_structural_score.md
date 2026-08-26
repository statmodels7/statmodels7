# The Score in a Structural Term's Own Parameters

The derivative of the weighted log-likelihood in the unconstrained
parameters of each structural term, one entry per term.

## Usage

``` r
statmod_structural_score(spec, coef, design = statmod_design(spec))
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design, whose structural state holds each term's current
  parameters.

## Value

A named list of numeric vectors, one entry per structural term, keyed by
the term's key in the specification, each as long as that term has
parameters. An empty list when the model carries no structural term.

## Details

The filter's Jacobian is already the **total** derivative of the
predictor in the term's own parameters, with the recursion propagated
inside it. So a chain rule over the observations is the whole of the
work here and no reverse pass is needed.

The reverse pass answers a different question: the derivative in the
coefficients of the equations, where a coefficient reaches the level
only through the scores at earlier times. That is
[`modelterms7::term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.html).

## See also

[`statmod_filter_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_filter_at.md)
for the filter run this reads,
[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
for the coefficients' own score.
