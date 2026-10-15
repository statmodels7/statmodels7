# Fit the Coefficients and a Filter's Parameters in One System

Minimizes the objective over the stacked coefficients of every equation
and the free parameters of a structural term of the filter shape in one
system, instead of alternating between the two blocks.

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

A list with `par` (the stacked coefficients only; the filter's free
parameters are written into the design's structural state), `value`,
`converged` and `iterations`.

## Details

The default optimizer is
[`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html),
with the exact gradient of both blocks and the exact observed
information over both together
([`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md),
plus the penalty's Hessian), run to its own stopping rule. An
alternation would hand the term's parameters to an optimizer of their
own at every pass, and every iteration of that optimizer would run the
recursion and its adjoint again with the coefficients held at a point
that is about to move.

The unknowns are ordered as the information orders them, the
coefficients first and the term's free parameters after, so no
permutation is needed; a level an intercept in the same equation carries
is held and leaves the system, exactly as it leaves the information.

The objective is evaluated once outside the guard, so that a term that
cannot be evaluated at all raises where it can be read; inside the
search a non-finite value belongs to the point (a filter whose loadings
put it outside the region where its recursion is bounded), and the
search steps back from that point instead of abandoning the run.

## See also

[`statmod_fit_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_structural.md),
which fits the term's parameters alone and is what a term of the
likelihood shape still uses.
