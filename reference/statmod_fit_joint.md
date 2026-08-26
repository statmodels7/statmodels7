# Fit the Coefficients and a Filter's Parameters in One System

One Newton step over the stacked coefficients of every equation and the
free parameters of a structural term of the filter shape, instead of
alternating between the two blocks.

## Usage

``` r
statmod_fit_joint(
  spec,
  design,
  obj,
  beta,
  hyper,
  optimizer = NULL,
  verbose = FALSE
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- obj:

  The objective, as
  [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  returns it.

- beta:

  The coefficients to start from.

- hyper:

  The hyperparameters.

- optimizer:

  An `optimizers7` optimizer, or `NULL` for Newton's method with the
  exact joint information.

- verbose:

  Logical; report the run.

## Value

A list with `par`, `value`, `converged` and `iterations`.

## Details

The alternation was not a statement about the model: the exact gradient
of both blocks and the exact observed information over both together are
already available, the second as
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md),
which was built for
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
and discarded for the fit. What the alternation cost is filter runs:
each sweep handed the term's parameters to an optimizer of their own,
and every one of that optimizer's iterations ran the recursion and its
adjoint again, with the coefficients held at a point that was about to
move.

The unknowns are ordered as the information orders them, the
coefficients first and the term's free parameters after, so no
permutation is needed; a level an intercept in the same equation carries
is held and leaves the system, exactly as it leaves the information.

The objective is evaluated once outside the guard, so that a term that
cannot be evaluated at all raises where it can be read; inside the
search a non-finite value is a statement about the point, a filter whose
loadings put it outside the region where its recursion is bounded, and
the search must step back from it, never abandon the run.

## See also

[`statmod_fit_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_structural.md),
which fits the term's parameters alone and is what a term of the
likelihood shape still uses.
