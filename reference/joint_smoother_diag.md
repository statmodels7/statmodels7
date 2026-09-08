# The Model's Smoother Over the Coefficients and a Filter's Parameters

The diagonal of \\F = (K+S)^{-1}K\\ on the vector a model carrying a
score-driven filter actually estimates: the coefficients of every
equation followed by the term's own free parameters.

## Usage

``` r
joint_smoother_diag(spec, coef, design, hyper)
```

## Arguments

- spec:

  The specification.

- coef:

  The coefficients.

- design:

  The design.

- hyper:

  The hyperparameters.

## Value

A list with `beta`, the diagonal over the stacked coefficients, `zeta`,
the diagonal over the structural term's free parameters, named as the
term names them, and `failed`. `NULL` where the model carries no joint
vector at all, and a list whose `failed` is `TRUE` where it carries one
that could not be read. The two are different answers: the first leaves
a caller free to fall back on a rule that is right for such a model, and
the second leaves it with no reading, so that the count is reported
missing rather than taken from a rule meaning something else.

## Details

The effective degrees of freedom are the trace of the model's smoother
and a term's share is the trace of its own diagonal block. Where every
unknown is a coefficient that matrix is \\(H+S)^{-1}H\\ and
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
reads it directly. A filter's parameters are estimated beside the
coefficients and may carry a penalty of their own, so for such a model
the same definition is read on the joint vector: \\K\\ is
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
and \\K+S\\ the matrix the criterion's own determinant is taken of,
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md).

## One apiece is this matrix's own answer, exactly

Writing \\P\\ for the coordinates some penalty covers, \\S\\ has a zero
row and a zero column outside \\P\\, so \$\$F = M^{-1}K = M^{-1}(M - S)
= I - M^{-1}S,\$\$ and for \\j\\ outside \\P\\ the whole column
\\S\_{\cdot j}\\ vanishes, which makes \\F\_{jj}\\ exactly one. The rule
this replaces – one degree of freedom per free parameter of a structural
term, on the reading that those are estimated and unpenalized – is
therefore what this matrix says wherever that rule is right, and
measured it is so to the last bit: the unpenalized coordinates come back
at 1.000000000000000 with a gap of 0.000e+00. What moves is a coordinate
a prior shrinks, which since 0.24.0 a filter's own parameters may be.
Measured on a developed loading over eight groups the eight deviations
read between 0.080 and 0.247 where the count says eight, and under a
covariance class over two developed parameters sixteen of them read
between -0.0003 and 0.829 where the count says sixteen.

A structural term that mixes over latent states rather than shifting the
predictor has no such joint determinant,
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
answering `NULL` for it, and falls back to the coefficient-only reading.

The joint route reads the OBSERVED information, there being no expected
one for a filter:
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
assembles from
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html),
and the criterion of such a model already reads that same matrix. A
model carrying no filter never reaches here and is untouched.

## See also

[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md),
its only caller.
