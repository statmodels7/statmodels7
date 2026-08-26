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

  The names of the coefficients `A` is indexed by.

## Value

The inverse.

## Details

A failure here is a statement about the fit and not about the
arithmetic: at a maximum the penalized information is positive definite,
so a matrix that is not says something about where the run stopped.

The test is `lmin > tol * ref` on the smallest eigenvalue, and never
whether [`chol()`](https://rdrr.io/r/base/chol.html) raised. On an
exactly singular matrix the pivot that should be zero comes out positive
or negative according to rounding, so the second answer differs between
platforms for one matrix.

`ref` is the matrix's own scale, or the scale of the unpenalized
information where the caller holds it. That is what separates a flat
direction from the scale separation a large smoothing parameter
legitimately produces. Returning a pseudo-inverse instead would give a
standard error for a direction the data does not identify.

The smallest eigenvalue is estimated, never computed, from LAPACK's
condition estimator (`dpocon`) read on the Cholesky factor the inverse
needs anyway: `rcond` is \\1/(\lVert A\rVert_1\lVert A^{-1}\rVert_1)\\,
so `rcond * ||A||_1` is \\1/\lVert A^{-1}\rVert_1\\, which for a
symmetric matrix lies between \\\lambda\_{\min}/\sqrt{p}\\ and
\\\lambda\_{\min}\\. The estimate therefore errs on the small side and
the test is conservative by at most a factor \\\sqrt{p}\\, plus whatever
the estimator's own slack is; the two cases it has to keep apart are
separated by some fifty orders of magnitude, so neither reaches the
other. It replaced a full eigendecomposition, which answers the same
question exactly and costs \\O(p^3)\\ with a large constant: measured at
\\p = 1022\\, 1.18 s against the Cholesky's 0.25.

The message names the directions, never the causes. A first version
offered two causes, the run not having reached a maximum and two columns
of the design carrying the same information, and on a Student t fitted
to `iris` neither was right: the design was full rank and the score was
4e-5. What had happened is the third and commonest case, a parameter
drifting to where its information vanishes, and no list of guesses would
have said so. The eigenvector of the smallest eigenvalue does: it is
read off and the coefficients that load on it are printed.
