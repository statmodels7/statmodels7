# Solve One Weighted Least Squares Step

Returns the increment \\\delta\\ solving \\(H + S)\delta = u\\ by the
requested decomposition.

## Usage

``` r
iwls_solve(pieces, u, how, damp = 0)
```

## Arguments

- pieces:

  A list with `R`, `C` and `A`, as
  [`iwls_pieces`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  builds it.

- u:

  The right-hand side.

- how:

  The decomposition.

- damp:

  The Levenberg damping \\\lambda\\, zero for the plain scoring step.

  **A coordinate at a boundary is held.** Where a parameter has reached
  the clamp its link keeps it strictly inside, the family's curvature
  there is zero or `NaN` and no step can move that coordinate. Such
  coordinates are dropped from the system and the rest is solved, which
  is the active set the boundary defines; the step stays a descent step
  for the reduced problem. Holding them is what keeps one boundary
  coordinate from stopping the others: with a single `NaN` in the
  curvature the whole step came back zero, and a Student t whose \\\nu\\
  had reached `double.xmax` stopped with a score of -617.6 in
  \\\sigma\\, an ordinary coordinate, while \\\nu\\'s own was exactly 0.

## Value

A list with `delta`, `rank`, `route` and `held`, the positions dropped
from the system.

## Details

`"qr"` and `"svd"` decompose the augmented matrix \\\[R;\\ C\]\\, whose
cross-product IS the penalized information, so \\X'X\\ is never formed
and the conditioning is never squared. That route needs the
per-observation curvature to be positive definite and the penalty
positive semidefinite; where either fails – an observed Hessian far from
the optimum, a non-convex penalty – there is no square root to take, and
the step falls back to the assembled matrix with its eigenvalues
floored, reporting in `route` that it did.

`"chol"` and `"chol_crossprod"` factor the assembled information
directly: faster per iteration and worse conditioned, which is the trade
the caller is choosing.

**The damping is Levenberg's and is zero unless a step has failed.**
`damp` adds \\\lambda I\\ to the system, which shortens a coordinate in
proportion to how little curvature it has: one with a diagonal of 0.24
beside neighbours at 2328 is shortened by \\(0.24+\lambda)/0.24\\ while
theirs move by \\(2328+\lambda)/2328\\, which is the differential shrink
a scalar step length cannot give. It is added to the augmented design as
further rows, so the QR route keeps its conditioning; the identity and
not \\\mathrm{diag}(K)\\, because a proportional damping shrinks every
coordinate alike and would leave the disparity where it was.
