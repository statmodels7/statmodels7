# The Information of the Weighted Log-Likelihood

Computes the negative Hessian in the coefficients, assembled block by
block from the distribution's own link-scale second derivatives.

## Usage

``` r
statmod_information_at(
  spec,
  coef,
  design = statmod_design(spec),
  expected = TRUE,
  emethod = NULL,
  H = NULL,
  index = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design, refreshed at `coef` if any term needs it.

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

  The second-derivative components at this point, as
  [`statmod_family_hessian()`](https://statmodels7.github.io/statmodels7/reference/statmod_family_hessian.md)
  returns them, or `NULL` to compute them here. A mixture over regimes
  reads a different set per component and ignores this.

- index:

  `NULL` for the whole matrix, or integer positions in the stacked
  coefficients: the submatrix over those positions, in the order given,
  assembled only from those columns.

## Value

A symmetric `p x p` matrix over the stacked coefficients, where `p` is
their total count across the equations. The matrix is a Matrix object
when any equation's design is sparse and a base matrix otherwise.

## Details

Block \\(a, b)\\ is \\X_a'\\\mathrm{diag}(w\\h\_{ab})\\X_b\\ with
\\h\_{ab}\\ the per-observation second derivative of the log-density in
\\(\eta_a, \eta_b)\\, negated. The blocks are assembled in the storage
type of the designs, so the information is sparse for a model with a
sparse equation.

`expected = TRUE` gives the expected information, which Fisher scoring
inverts and which is positive definite wherever the family is regular.
`expected = FALSE` gives the negated observed Hessian, which far from
the optimum may be indefinite.

`emethod` is the `method` of
[`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html)
and is read only when `expected` is `TRUE`.

With `index` only the requested columns of each design enter a cross
product, so that reading a few rows and columns does not require
\\X'WX\\ over every coefficient.

For a likelihood mixed over latent states the matrix is the
complete-data information, the ordinary information averaged over the
smoothed states. It is positive definite and serves as a scoring matrix,
and it is not the observed information of the mixture. The variance of a
fit comes from
[`statmod_regime_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_regime_information.md).

## See also

[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
for the first derivative,
[`info_blocks()`](https://statmodels7.github.io/statmodels7/reference/info_blocks.md)
for the per-observation blocks that a square-root route uses instead,
and
[`statmod_family_hessian()`](https://statmodels7.github.io/statmodels7/reference/statmod_family_hessian.md)
for the components both read.
