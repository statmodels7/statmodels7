# Solve One Weighted Least Squares Step

Returns the increment \\\delta\\ solving \\(H + S)\delta = u\\ by the
requested decomposition.

## Usage

``` r
iwls_solve(pieces, u, how, damp = 0, frozen = integer(0))
```

## Arguments

- pieces:

  A list with `R`, `C` and `A`, as
  [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  builds it. Which of the three are present depends on `how`.

- u:

  The right-hand side \\X'g - S\beta\\, a numeric vector of length `p`.

- how:

  The decomposition: `"qr"`, `"svd"` or `"chol"`.

- damp:

  The Levenberg damping \\\lambda\\, a single non-negative number. `0`,
  the default, is the plain scoring step.

- frozen:

  Positions that the caller holds at their current values, an integer
  vector. They are dropped from the system exactly as a coordinate whose
  curvature is not finite is, so the step is that of the problem reduced
  to the rest. The default is empty.

## Value

A list of five:

- `delta`:

  the increment, a numeric vector of length `p`, zero in every held
  coordinate.

- `rank`:

  the rank of the system solved, an integer.

- `route`:

  which decomposition was used: the one requested, or `"qr->chol"` or
  `"svd->chol"` when the augmented route fell back to the assembled
  matrix with its eigenvalues floored.

- `held`:

  the positions dropped from the system, an integer vector, empty when
  every coordinate was free.

- `dropped`:

  the positions that the pivoted
  [`qr()`](https://rdrr.io/r/base/qr.html) left out as aliased, an
  integer vector, empty on every other route.

## The two families of route

`"qr"` and `"svd"` decompose the augmented matrix \\\[R;\\ C\]\\, whose
cross-product is the penalized information, so \\X'X\\ is never formed
and the condition number is never squared. That route needs the
per-observation curvature to be positive definite and the penalty
positive semidefinite. Where either fails, as with an observed Hessian
far from the optimum or a non-convex penalty, there is no square root to
take: the step falls back to the assembled matrix with its eigenvalues
floored, and records the fallback in `route` as `"qr->chol"` or
`"svd->chol"`.

`"chol"` factors the assembled information directly, after flooring the
eigenvalues of a matrix that is not positive definite.

## Levenberg damping

`damp` adds \\\lambda I\\ to the system, which shortens a coordinate in
proportion to how little curvature it has. A coordinate whose diagonal
is \\d\\ is shortened by \\(d + \lambda)/d\\, so that a coordinate with
a small diagonal is shortened more than one with a large diagonal: a
differential shrink, which a scalar step length cannot give. The damping
is zero unless a step has failed.

On the augmented route the damping enters as further rows of the
augmented design, so the QR route keeps its conditioning; on the
assembled route it is added to the diagonal. The matrix added is the
identity. A damping proportional to \\\mathrm{diag}(K)\\ would shrink
every coordinate alike and leave the disparity exactly where it was.

## A coordinate at a boundary is held

Where a parameter has reached the clamp with which its link keeps it
strictly inside its bounds, the family's curvature there is zero or
`NaN`, and that coordinate cannot be moved by any step. Those
coordinates are dropped from the system and the rest is solved, which is
the active set that the boundary defines; the step stays a descent step
for the reduced problem. Coordinates in `frozen` and coordinates with a
non-finite right-hand side are held in the same way.

Holding them keeps a single boundary coordinate from preventing the
others from moving: a `NaN` in one entry of the curvature would
otherwise make the whole step zero.

On the augmented route a coordinate is held when its column of the
square-root design or of the penalty factor has a non-finite entry. On
the assembled route the test is on the diagonal of the curvature. A
boundary coordinate makes its entire row non-finite, cross terms
included, so a test over whole columns would mark its neighbors too and
hold coordinates that can move. After the diagonal test, the columns of
the reduced matrix are tested again, and the test repeats until every
column is finite.

## See also

[`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
for the input,
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the `decomposition` argument that this function serves.
