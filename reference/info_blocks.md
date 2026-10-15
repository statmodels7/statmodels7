# The Per-Observation Information Blocks

Assembles \\\Omega_i\\, the \\K \times K\\ information of observation
\\i\\ with respect to the \\K\\ link-scale predictors, multiplied by
that observation's prior weight. This is the per-observation curvature a
scoring step is built from, before any design enters.

## Usage

``` r
info_blocks(spec, theta, expected = TRUE, approx = "opg", H = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  read for the distribution, the weights and the thread count.

- theta:

  The per-observation parameters on the parameter scale, a named list of
  \\n\\-vectors, one per distribution parameter.

- expected:

  `TRUE` for the expected information, `FALSE` for the negated observed
  Hessian.

- approx:

  How the expected information is approximated for a family that has no
  closed form, passed through to distributions7: `"opg"` (the default),
  `"bartlett"`, `"integrate"` or `"mc"`. Ignored when `expected` is
  `FALSE` or when the family computes its expected information exactly.

- H:

  The family's link-scale second derivatives, as
  [`statmod_family_hessian()`](https://statmodels7.github.io/statmodels7/reference/statmod_family_hessian.md)
  returns them, when the caller has already computed them at `theta`;
  `NULL` computes them here.

## Value

An \\n \times K \times K\\ numeric array, symmetric in its last two
indices, with \\K\\ the number of distribution parameters.

## Details

The derivatives come from the family on the link scale, so the chain
rule onto \\\eta\\ has already been applied by distributions7 and
nothing here multiplies by a link's derivative. With `expected = TRUE`
the expected information is taken and the blocks are positive definite
wherever the family is regular; with `expected = FALSE` the observed
Hessian is negated, which far from the optimum is routinely indefinite.
In that case
[`chol_blocks()`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md)
returns `NULL`.

The weights enter as given, without normalization, so a weight of two
counts an observation twice.

## See also

[`chol_blocks()`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md),
which factorizes these,
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
for the assembled \\Z'\Omega Z\\.
