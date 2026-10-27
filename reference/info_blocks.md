# The Per-Observation Information Blocks

Assembles \\\Omega_i\\, the \\K \times K\\ information of observation
\\i\\ with respect to the \\K\\ link-scale predictors, multiplied by
that observation's prior weight. This is the per-observation curvature
from which a scoring step is built, before any design enters.

## Usage

``` r
info_blocks(spec, theta, expected = TRUE, emethod = NULL, H = NULL)
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

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

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
rule onto \\\eta\\ has already been applied by distributions7 and the
function does not multiply by a link's derivative. With
`expected = TRUE` the expected information is taken and the blocks are
positive definite wherever the family is regular; with
`expected = FALSE` the observed Hessian is negated, which far from the
optimum is routinely indefinite. In that case
[`chol_blocks()`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md)
returns `NULL`.

The weights enter as given, without normalization, so a weight of two
counts an observation twice.

## See also

[`chol_blocks()`](https://statmodels7.github.io/statmodels7/reference/chol_blocks.md),
which factorizes these,
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
for the assembled \\Z'\Omega Z\\.
