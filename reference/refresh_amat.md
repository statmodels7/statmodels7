# The Weight for Contracting a Refreshable Block

Computes \\A_a\[i,j\] = w_i\sum_b c\_{ab}(i)\\(X_b
M\_{ab}^\top)\[i,j\]\\, the matrix with which a term's own derivative is
paired.

## Usage

``` r
refresh_amat(spec, design, M, params, npar, offs, ra, cw)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- M:

  The matrix against which the trace is taken.

- params, npar, offs:

  The block bookkeeping.

- ra:

  The term's rows in the stacked coefficient vector.

- cw:

  One length-\\n\\ vector per distribution parameter, or `NULL` where
  that parameter carries no coefficient.

## Value

A numeric matrix, one row per observation and one column per coefficient
of the term.

## Details

The same shape serves three quantities and differs only in \\c\_{ab}\\:
the curvature \\\ell\_{ab}\\ for the gradient's contraction, and the
third derivative already contracted in a direction for each of the two
mixed terms of the Hessian.
