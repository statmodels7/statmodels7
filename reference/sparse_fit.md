# Fit One Non-Smooth Block While the Others Are Held Fixed

Estimates the coefficients of one block with a kinked penalty while the
other coefficients are held fixed. A penalty over a structural term's
own parameters is fitted by
[`sparse_fit_structural()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit_structural.md).
Otherwise a compiled coordinate descent,
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md),
is used where the penalty describes its proximal operator as a table,
and an accelerated proximal gradient iteration runs on the block where
it does not, in which the smooth part of the objective supplies the
gradient and the penalty supplies its proximal operator.

## Usage

``` r
sparse_fit(
  obj,
  beta,
  block,
  hyper,
  maxit = 500,
  tol = 1e-08,
  verbose = FALSE,
  spec = NULL,
  design = NULL,
  expected = TRUE,
  approx = "opg"
)
```

## Arguments

- obj:

  The full objective.

- beta:

  The current stacked coefficients.

- block:

  One entry of `statmod_blocks()$sparse`.

- hyper:

  The hyperparameters.

- maxit:

  The iteration budget of the proximal gradient route and of the
  structural route. The coordinate descent route has a budget of its
  own, set in
  [`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
  and
  [`coord_call()`](https://statmodels7.github.io/statmodels7/reference/coord_call.md).

- tol:

  The stopping tolerance.

- verbose:

  Whether the optimizer of the proximal routes prints its own trace. The
  coordinate descent route does not use it.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).
  If `spec` or `design` is `NULL`, the coordinate descent is skipped.

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- expected:

  `TRUE` for the expected information in the working weights of the
  coordinate descent, `FALSE` for the observed one.

- approx:

  How the expected information is approximated for a family with no
  closed form, as in
  [`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md).

## Value

A list with `par` (the whole vector, updated in this block), `value`,
`converged`, `iterations` and `method` (`"coordinate descent"`,
`"proximal gradient"` or `"proximal gradient (structural)"`).

## Details

The coordinate descent of
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
reads the columns of the design and the running residual, which belong
to the model, so it is called here and not through an optimizer's
black-box interface. Where
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
returns `NULL` (the penalty has no proximal table at the steps of this
block, or the working weights are unusable), or where `spec` or `design`
is `NULL`, the route is
[`optimizers7::prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.html)
with
[`penalties7::penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.html),
accelerated.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
