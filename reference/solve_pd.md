# Invert a Matrix That Ought to Be Positive Definite

An inverse through the Cholesky factor, signaling an error naming the
matrix when a direction is flat.

## Usage

``` r
solve_pd(A, what, labels = NULL)
```

## Arguments

- A:

  A square matrix.

- what:

  What the matrix is, for the message.

- labels:

  The names of the coefficients that index `A`.

## Value

The inverse.

## Details

A failure here is a statement about the fit and not about the
arithmetic: at a maximum the penalized information is positive definite,
so a matrix that is not positive definite indicates that the run stopped
away from a maximum.

The matrix is first scaled to a unit diagonal, \\A_e = D^{-1/2} A
D^{-1/2}\\ with \\D\\ the diagonal of \\A\\. The test is
`lmin > 1e-12 * ||A_e||_1` on the smallest eigenvalue of \\A_e\\, and
never whether [`chol()`](https://rdrr.io/r/base/chol.html) raised. On an
exactly singular matrix the pivot that should be zero comes out positive
or negative according to rounding, so the second answer differs between
platforms for one matrix.

The scaling separates a flat direction from the scale separation that a
large smoothing parameter legitimately produces, because a direction is
flat only against its own scale. Returning a pseudo-inverse instead
would give a standard error for a direction that the data do not
identify.

The smallest eigenvalue is estimated and not computed, from LAPACK's
condition estimator (`dpocon`) read on the Cholesky factor that the
inverse needs anyway: `rcond` is \\1/(\lVert A\rVert_1\lVert
A^{-1}\rVert_1)\\, so `rcond * ||A||_1` is \\1/\lVert A^{-1}\rVert_1\\.
For a symmetric matrix the exact value of this product lies between
\\\lambda\_{\min}/\sqrt{p}\\ and \\\lambda\_{\min}\\. The estimator of
\\\lVert A^{-1}\rVert_1\\ is a lower bound, so the value read can exceed
the exact product, and can exceed \\\lambda\_{\min}\\ by a small factor.
The two cases to be separated differ by many orders of magnitude, so
neither reaches the other. The estimate replaces a full
eigendecomposition, which answers the same question exactly and costs
\\O(p^3)\\ with a larger constant than the Cholesky factorization.

The message names the coefficients that load on the eigenvector of the
smallest eigenvalue, and no cause. A flat direction can arise from a run
that has not reached a maximum, from two columns of the design carrying
the same information, or from a parameter that has drifted to where its
information vanishes.
