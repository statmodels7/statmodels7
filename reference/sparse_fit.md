# Fit One Non-Smooth Block, the Others Held Fixed

Runs a proximal gradient iteration on the block's own coefficients, the
smooth part of the objective supplying the gradient and the penalty its
proximal operator.

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

  The iteration budget.

- tol:

  The stopping tolerance.

- verbose:

  Whether the optimizer prints its own trace.

## Value

A list with `par` (the whole vector, updated in this block), `value`,
`converged` and `iterations`.

## Details

The route is
[`optimizers7::prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.html)
with
[`penalties7::penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.html),
accelerated: measured, at a condition number of 3 the plain iteration
wins narrowly, at 55 it is 4153 iterations against 126, and at 480 the
plain one does not converge in 50000. A coordinate descent on the same
objective is 1.1 to 5.3 times faster again and is the next thing to
write. It needs the columns of the design and the running residual,
which is the model itself, so it belongs here and could not live behind
an optimizer's black-box interface.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
