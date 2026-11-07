# The Compiled Alternation of a Model with a Kinked Penalty

`engine_alternate_fit()` runs
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
in compiled code (`src/engine_kink.h`) for a model whose penalty has a
kink: at every pass the curvature of a scaled SCAD or MCP, the smooth
block fitted on a subset as
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
fits it, every kinked block fitted by the coordinate descent of
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md),
and at the end of the passes the reading of
[`alternation_settled()`](https://statmodels7.github.io/statmodels7/reference/alternation_readings.md).
It returns the list that
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
returns, or `NULL`, in which case
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
runs the passes in R from the same coefficients.

## Usage

``` r
engine_alternate_fit(
  spec,
  design,
  blocks,
  hyper,
  inner_optimizer,
  beta,
  expected,
  emethod,
  maxit,
  tol,
  vb,
  obj
)

engine_alternate_blocks(spec, design, blocks, units, obj, beta, hf)
```

## Arguments

- spec, design, blocks, hyper, inner_optimizer, beta, expected, emethod,
  maxit, tol, vb:

  As for
  [`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md).

- obj:

  The objective at the start,
  [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md).

- units:

  The penalized units,
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md).

- hf:

  The held positions,
  [`held_positions()`](https://statmodels7.github.io/statmodels7/reference/held_positions.md).

## Value

`engine_alternate_fit()`: the list of
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md),
or `NULL`. `engine_alternate_blocks()`: a list of block descriptions, or
`NULL` where a block is not covered.

## Details

A model is covered when
[`engine_route()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md)
covers its smooth block with the closed-form expected information and
`emethod` is `NULL`; nothing is printed (`vb$blocks` and `vb$inner` are
`FALSE`); the model has one thread, a dense design and no structural,
moving or frozen term; and every kinked block is a lasso, an elastic
net, a SCAD or an MCP over a dense block, with a diagonal map or none,
whose coordinates no other penalized unit covers and none of which is
held. Where the compiled alternation meets a branch that it does not
reproduce (the proximal fallback of
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md))
or a point where the R route signals an error, it returns `NULL`.

`engine_alternate_blocks()` describes the kinked blocks to the compiled
alternation: for each, its penalized unit, its equation, its columns and
the other columns of the equation, whether it is solved with the
equation's intercept profiled out
([`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md)),
whether
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
writes its curvature and with which centering, the previous point's kink
on a path, and the columns themselves.

## See also

[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md),
[`engine_inner_fit()`](https://statmodels7.github.io/statmodels7/reference/engine_inner_fit.md)
