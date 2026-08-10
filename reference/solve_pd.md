# Invert a Matrix That Ought to Be Positive Definite

A Cholesky inverse, signalling an error naming the matrix when the
factor does not exist.

## Usage

``` r
solve_pd(A, what)
```

## Arguments

- A:

  A square matrix.

- what:

  What the matrix is, for the message.

## Value

The inverse.

## Details

A failure here is a statement about the fit rather than about the
arithmetic: at a maximum the penalized information is positive definite,
so a factor that does not exist says the run stopped somewhere that is
not one, or that two columns of the design carry the same information.
Returning a pseudo-inverse instead would give a standard error for a
direction the data does not identify.
