# The Observed Information Over the Coefficients and a Filter's Parameters

The information of a model carrying a structural filter, over the
coefficients of every equation and the term's own parameters, with the
recursion's own curvature in it.

## Usage

``` r
statmod_full_information(spec, coef, design = statmod_design(spec))

statmod_full_information_impl(spec, coef, design, params, ev)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

- params, ev:

  The parameter names and the evaluated predictors, already in hand at
  the caller.

## Value

A square matrix over the coefficients and the term's parameters, or the
ordinary information where there is no structural term.

## Details

The gradient of such a model is exact already, through
[`modelterms7::term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.html).
Its curvature is not: the matrix the scoring step inverts is assembled
as though the level were an offset, which is a legitimate scoring matrix
and is not the information. Writing \\u\\ for the coefficients followed
by the term's parameters on the unconstrained scale, \\V\_{q,t}\\ for
the derivative of equation \\q\\'s predictor in \\u\\, and \\E_t\\ for
the second derivative of the predictor the filter produces,

\$\$-\frac{\partial^2 \ell}{\partial u \partial u^\top} = -\sum_t w_t
\sum\_{q,r} \ell\_{qr,t} V\_{q,t}^\top V\_{r,t} - \sum_t w_t \ell\_{p,t}
E_t.\$\$

Only the equation carrying the filter has a \\V\\ that is not its own
design: there it is the forward Jacobian of the recursion, which
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
returns beside the contracted \\E\\. The third derivatives the second
sum needs are distributions7's, in closed form for every family.

The term's parameters are the last columns, which is the convention
`term_curvature()` shares with its caller.

A term that mixes over latent states, in place of shifting the
predictor, is routed to
[`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md),
whose Hessian comes from the same forward recursion the likelihood does;
a model carrying neither gets the ordinary information.

## See also

[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md),
[`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md)
