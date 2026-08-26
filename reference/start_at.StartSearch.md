# Starting Values From a Global Search

Runs the optimizer once on the likelihood, over the non-convex
coordinates and, where the model has one, over a structural term's own
parameters.

## Arguments

- strategy:

  A `StartSearch` object.

- spec, design, obj, ...:

  As in
  [`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md).

## Value

A named list of numeric vectors.

## Details

A structural term's parameters do not live in the coefficient vector;
they are held in the design's structural state, and the objective reads
them from there. The search therefore sets them into that state and
evaluates the likelihood, exactly as
[`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
does, and leaves the best it found in place — which is how a strategy
returning coefficients can nonetheless start a filter somewhere better
than its own `term_start()`.
