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
  [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  builds it. Which of the three are present depends on `how`.

- u:

  The right-hand side \\X'g - S\beta\\, a numeric vector of length `p`.

- how:

  The decomposition: `"qr"`, `"svd"`, `"chol"` or `"chol_crossprod"`.

- damp:

  The Levenberg damping \\\lambda\\, a single non-negative number. `0`,
  the default, is the plain scoring step.

## Value

A list of four:

- `delta`:

  the increment, a numeric vector of length `p`, zero in every held
  coordinate.

- `rank`:

  the rank of the system solved, an integer.

- `route`:

  which decomposition was used, including `"floor"` when the augmented
  route fell back to the floored assembled matrix.

- `held`:

  the positions dropped from the system, an integer vector, empty when
  nothing was held.

## The two families of route

`"qr"` and `"svd"` decompose the augmented matrix \\\[R;\\ C\]\\, whose
cross-product is the penalized information, so \\X'X\\ is never formed
and the condition number is never squared. That route needs the
per-observation curvature to be positive definite and the penalty
positive semidefinite. Where either fails, as with an observed Hessian
far from the optimum or a non-convex penalty, there is no square root to
take: the step falls back to the assembled matrix with its eigenvalues
floored, and says so in `route`.

`"chol"` and `"chol_crossprod"` factor the assembled information
directly. Faster per iteration, worse conditioned, and the trade the
caller chose.

## The damping is Levenberg's, and is zero unless a step has failed

`damp` adds \\\lambda I\\ to the system, which shortens a coordinate in
proportion to how little curvature it has. A coordinate whose diagonal
is 0.24 beside neighbors at 2328 is shortened by \\(0.24 +
\lambda)/0.24\\ while theirs move by \\(2328 + \lambda)/2328\\: a
differential shrink, which a scalar step length cannot give.

It goes in as further rows of the augmented design, so the QR route
keeps its conditioning. The matrix added is the identity. A damping
proportional to \\\mathrm{diag}(K)\\ would shrink every coordinate alike
and leave the disparity exactly where it was.

## A coordinate at a boundary is held

Where a parameter has reached the clamp its link keeps it strictly
inside, the family's curvature there is zero or `NaN`, and no step can
move that coordinate. Those coordinates are dropped from the system and
the rest is solved, which is the active set the boundary defines; the
step stays a descent step for the reduced problem.

Holding them keeps one boundary coordinate from stopping the others.
With a single `NaN` in the curvature the whole step came back zero, and
a Student t whose \\\nu\\ had reached `double.xmax` stopped with a score
of -617.6 in \\\sigma\\, an ordinary coordinate nowhere near stationary,
while \\\nu\\'s own score was exactly 0.

The test is on the **diagonal** of the curvature. A boundary coordinate
makes its entire row non-finite, cross terms included, so a test over
whole columns marks its neighbors too and holds coordinates that were
perfectly movable.

## See also

[`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
for the input,
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the `decomposition` argument this serves.
