# Fit the Smooth Block by Iterated Weighted Least Squares

Runs the scoring iteration on the objective of
[`statmod_objective`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
at fixed hyperparameters, returning the coefficients and the history of
the run.

## Usage

``` r
iwls_fit(obj, start, method, n, pieces_at, verbose = FALSE)
```

## Arguments

- obj:

  The objective, from
  [`statmod_objective`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md).

- start:

  The starting coefficients, stacked.

- method:

  An
  [`Iwls`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

- n:

  The number of observations, for the scaled stopping rule.

- pieces_at:

  A function of the stacked coefficients returning what
  [`iwls_solve`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
  needs, as
  [`iwls_pieces`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  builds it.

- verbose:

  Whether to print a line per iteration.

## Value

A list with `par`, `value`, `converged`, `iterations`, `score` and
`history`.

## Details

Two things the loop does that a plain Newton iteration does not, both
recorded elsewhere in the toolkit as the reason a run reported failure
at the answer. A non-positive-definite curvature is repaired by flooring
its eigenvalues rather than abandoning the start, since
[`solve()`](https://rdrr.io/r/base/solve.html) would otherwise force
one. And the stopping rule is read at the ITERATE, on a score scaled by
the sample size, so that a threshold means the same thing at \\n = 10\\
and at \\n = 10^7\\ while the objective itself stays unaveraged.

## See also

[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
