# The Model's Smoother Over the Coefficients and a Filter's Parameters

Returns the diagonal of \\F = (K+S)^{-1}K\\ on the vector that a model
carrying a score-driven filter actually estimates: the coefficients of
every equation followed by the term's own free parameters.

## Usage

``` r
joint_smoother_diag(
  spec,
  coef,
  design,
  hyper,
  expected = FALSE,
  emethod = NULL
)
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

- expected, emethod:

  Whether the matrix is the expected information conditional on the
  past, and its `method`.

## Value

A list with `beta`, the diagonal over the stacked coefficients, `zeta`,
the diagonal over the structural term's free parameters, named as the
term names them, and `failed`. The result is `NULL` where the model
carries no joint vector at all, and a list whose `failed` is `TRUE`
where it carries one that could not be read. The two results differ:
with the first, the caller may fall back on a rule that is right for
such a model, and with the second, the caller has no reading, so that
the count is reported missing rather than taken from a rule that means
something else.

## Details

The effective degrees of freedom are the trace of the model's smoother
and a term's share is the trace of its own diagonal block. Where every
unknown is a coefficient that matrix is \\(H+S)^{-1}H\\ and
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
reads it directly. A filter's parameters are estimated beside the
coefficients and may carry a penalty of their own, so for such a model
the same definition is read on the joint vector: \\K\\ is
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
and \\K+S\\ the matrix whose determinant the criterion takes,
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md).

## The count of an unpenalized parameter

Writing \\P\\ for the coordinates that some penalty covers, \\S\\ has a
zero row and a zero column outside \\P\\, so \$\$F = M^{-1}K =
M^{-1}(M - S) = I - M^{-1}S,\$\$ and for \\j\\ outside \\P\\ the whole
column \\S\_{\cdot j}\\ vanishes, which makes \\F\_{jj}\\ exactly one.
The rule that counts one degree of freedom per free parameter of a
structural term, on the reading that those parameters are estimated and
unpenalized, therefore agrees with this matrix wherever that rule is
right. The count differs for a coordinate that a prior shrinks, which a
filter's own parameters may be (for example the deviations of a
developed loading, or the parameters under a covariance class); such a
coordinate counts less than one.

A structural term that mixes over latent states rather than shifting the
predictor has no such joint determinant, because
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
returns `NULL` for it, so the count falls back to the coefficient-only
reading.

The joint route reads the observed information by default, and the
expected information conditional on the past with `expected = TRUE` (see
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md));
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
assembles from
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html),
and the criterion of such a model already reads that same matrix. A
model carrying no filter never reaches here and is untouched.

## See also

[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
and
[`summary_class_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_class_blocks.md),
its callers.
