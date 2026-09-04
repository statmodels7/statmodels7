# The Information of the Weighted Log-Likelihood

The negative Hessian in the coefficients, assembled block by block from
the distribution's own link-scale second derivatives.

## Usage

``` r
statmod_information_at(
  spec,
  coef,
  design = statmod_design(spec),
  expected = TRUE,
  approx = "opg",
  H = NULL
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

- approx:

  How the expected information is approximated for a family with no
  closed form: `"bartlett"`, `"integrate"` or `"mc"`.

- H:

  The second-derivative components at this point, as
  [`statmod_family_hessian()`](https://statmodels7.github.io/statmodels7/reference/statmod_family_hessian.md)
  returns them, or `NULL` to ask for them. A mixture over regimes reads
  a different set per component and ignores this.

## Value

A symmetric `p x p` matrix over the stacked coefficients, `p` being
their total count across the equations. A Matrix object when any
equation's design is sparse, a base matrix otherwise.

## Details

Block \\(a, b)\\ is \\X_a'\\\mathrm{diag}(w\\h\_{ab})\\X_b\\ with
\\h\_{ab}\\ the per-observation second derivative of the log-density in
\\(\eta_a, \eta_b)\\, negated. The blocks are assembled in the storage
the designs call for, so a model with a sparse equation gets a sparse
information.

`expected = TRUE` gives the expected information, which Fisher scoring
inverts and which is positive definite wherever the family is regular.
`expected = FALSE` gives the negated observed Hessian, which far from
the optimum may be indefinite.

`approx` reaches distributions7 and is read only where the family has no
closed expected information.

## See also

[`statmod_score_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_at.md)
for the first derivative,
[`info_blocks()`](https://statmodels7.github.io/statmodels7/reference/info_blocks.md)
for the per-observation blocks a square-root route uses instead, and
[`statmod_family_hessian()`](https://statmodels7.github.io/statmodels7/reference/statmod_family_hessian.md)
for the components both read.
