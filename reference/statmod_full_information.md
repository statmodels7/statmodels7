# The Information Over the Coefficients and a Filter's Parameters

Computes the information of a model carrying a structural filter, over
the coefficients of every equation and the term's own parameters: the
observed information, with the recursion's own curvature in it, or the
expected information conditional on the past.

## Usage

``` r
statmod_full_information(
  spec,
  coef,
  design = statmod_design(spec),
  expected = FALSE,
  emethod = NULL
)

statmod_full_information_impl(
  spec,
  coef,
  design,
  params,
  ev,
  expected = FALSE,
  emethod = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

- expected:

  `FALSE` (the default) for the observed information, `TRUE` for the
  expected information conditional on the past.

- emethod:

  The `method` of
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html),
  read only where `expected` is `TRUE`.

- params, ev:

  The parameter names and the evaluated predictors, already in hand at
  the caller.

## Value

A square matrix over the coefficients followed by the term's free
parameters (a level held because an intercept carries it leaves the
matrix), or the ordinary information where there is no structural term.

## Details

The gradient of such a model is exact already, through
[`modelterms7::term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.html).
The matrix that the scoring step inverts is assembled as though the
level were an offset; it serves as a scoring matrix and differs from the
observed information, which this function computes. Writing \\u\\ for
the coefficients followed by the term's parameters on the unconstrained
scale, \\V\_{q,t}\\ for the derivative of equation \\q\\'s predictor in
\\u\\, and \\E_t\\ for the second derivative of the predictor that the
filter produces,

\$\$-\frac{\partial^2 \ell}{\partial u \partial u^\top} = -\sum_t w_t
\sum\_{q,r} \ell\_{qr,t} V\_{q,t}^\top V\_{r,t} - \sum_t w_t \ell\_{p,t}
E_t.\$\$

Only the equation carrying the filter has a \\V\\ that is not its own
design: there it is the forward Jacobian of the recursion, which
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
returns beside the contracted \\E\\. The third derivatives that the
second sum needs are those of distributions7, from
[`distributions7::distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.html).

The term's parameters are the last columns, which is the convention that
`term_curvature()` shares with its caller.

With `expected = TRUE` the second sum is dropped and \\\ell\_{qr,t}\\ is
replaced by its expectation given the observations before \\t\\:

\$\$\mathcal{I}\_{\rm cond} = \sum_t w_t \sum\_{q,r}
\mathcal{I}\_{qr}(\eta_t)\\ V\_{q,t}^\top V\_{r,t},\$\$

where \\\mathcal{I}\_{qr}(\eta_t)\\ is the expected information of the
family on the link scale at the predictors of observation \\t\\, from
[`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html)
with `method = emethod`. The predictors \\\eta_t\\ and their Jacobian
\\V\_{q,t}\\ depend on the observations before \\t\\ only, and the score
at \\t\\ has conditional mean zero, so the second sum has conditional
mean zero and drops out. The matrix is positive semidefinite and needs
the first derivative of the recursion only. With
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
\\\mathcal{I}\_{qr}(\eta_t)\\ is the product of the scores at the data,
and the matrix is the sum of the outer products of the per-observation
scores in \\u\\.

A term that mixes over latent states, in place of shifting the
predictor, is routed to
[`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md),
whose Hessian comes from the same forward recursion that computes the
likelihood; for a model that carries neither kind of term, the ordinary
information is used. The likelihood of a model with a term of that kind
is a mixture over latent components, and the expectation of its second
derivatives has no closed form, so `expected = TRUE` signals an error
there (see
[`structural_expected_check()`](https://statmodels7.github.io/statmodels7/reference/structural_expected_check.md)).

## See also

[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md),
[`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md)
