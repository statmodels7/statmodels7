# The Two Readings That a Summary Prints About a Point

Computes the largest gradient in standard-error units, \\\max_j
\|g_j\|/\sqrt{\|A\_{jj}\|}\\, and the smallest eigenvalue of the
equilibrated curvature \\D^{-1/2} A D^{-1/2}\\ with \\D =
\|\mathrm{diag}(A)\|\\.

## Usage

``` r
certificate_readings(g, A, keep = seq_along(g))
```

## Arguments

- g:

  The gradient.

- A:

  The curvature, a matrix of the same order.

- keep:

  The coordinates to read.

## Value

A named numeric vector, `gradient` and `eigen`, `NA` where there is
nothing to read.

## Details

Both are dimensionless, so both are independent of the units of a
coefficient and of the sample size. \\\|g_j\|/\sqrt{A\_{jj}}\\ is the
length of the one-coordinate Newton step measured in that coordinate's
standard error conditional on the others, \\1/\sqrt{A\_{jj}}\\, which is
smaller than the marginal one; equivalently it is \\\sqrt{2\Delta_j}\\,
with \\\Delta_j = g_j^2/(2A\_{jj})\\ the amount that this one-coordinate
step would still gain in the objective's own units. The equilibrated
matrix has a unit diagonal, so its smallest eigenvalue is at most one:
near zero it reports a direction that the curvature barely identifies,
and below zero a point that is not an optimum.

`A` is the curvature with the sign that makes an optimum positive
definite: the penalized information for the inner fit, the negated outer
Hessian for the criterion. A coordinate with zero curvature reads an
infinite gradient where its gradient is not zero, and keeps a unit scale
in the equilibration.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md)
