# Invert a Matrix That Ought to Be Positive Definite

A Cholesky inverse, signalling an error naming the matrix when the
factor does not exist.

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

A failure here is a statement about the fit rather than about the
arithmetic: at a maximum the penalized information is positive definite,
so a matrix that is not says something about where the run stopped. The
test is `min(ev) > tol * max(ev)` on the eigenvalues rather than whether
[`chol()`](https://rdrr.io/r/base/chol.html) raised, because on an
exactly singular matrix the latter is decided by rounding and differs
between platforms. Returning a pseudo-inverse instead would give a
standard error for a direction the data does not identify.

The message names the directions rather than the causes. A first version
offered two – the run not having reached a maximum, or two columns of
the design carrying the same information – and on a Student t fitted to
iris NEITHER was right: the design was full rank and the score was 4e-5.
What had happened is the third and commonest case, a parameter drifting
to where its information vanishes, and no list of guesses would have
said so. The eigenvector of the smallest eigenvalue does: it is read off
and the coefficients that load on it are printed.
