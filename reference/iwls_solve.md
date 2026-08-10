# Solve One Weighted Least Squares Step

Returns the increment \\\delta\\ solving \\(H + S)\delta = u\\ by the
requested decomposition.

## Usage

``` r
iwls_solve(pieces, u, how)
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

## Value

A list with `delta`, `rank` and `route`.

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
