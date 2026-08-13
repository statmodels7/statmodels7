# The Alternation Between the Smooth Block and the Rest

Fits the terms whose penalties are twice differentiable in one system
and each remaining block by a method of its own, sweeping until the
objective stops moving.

## Usage

``` r
statmod_alternate(
  spec,
  design,
  blocks,
  hyper,
  inner_optimizer,
  beta,
  expected,
  approx,
  maxit,
  tol,
  vb
)
```

## Arguments

- spec:

  The specification.

- design:

  The design.

- blocks:

  The block split.

- hyper:

  The hyperparameters.

- inner_optimizer:

  How the smooth block is fitted.

- beta:

  The starting coefficients, stacked.

- expected:

  Whether the information is the expected one.

- approx:

  The approximation for the expected information.

- maxit, tol:

  The budget and the tolerance.

- vb:

  The resolved verbosity.

## Value

A list with `par`, `value`, `converged`, `obj`, `hist_blocks` and
`hist_inner`.

## Details

It is a function of its own because the outer search calls it once per
hyperparameter it tries, warm-started from the previous coefficients.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
