# The Flat Directions of the Outer Curvature

The directions along which the equilibrated outer curvature \\A/(ss')\\,
\\s_j = \sqrt{\|A\_{jj}\|}\\, has an eigenvalue of at most `flat` in
absolute value, over the coordinates `keep`.

## Usage

``` r
certificate_ridges(g, A, keep, flat, tol)
```

## Arguments

- g:

  The outer gradient.

- A:

  The outer curvature, positive at a maximum.

- keep:

  The coordinates to read, those not already reported at an edge.

- flat:

  The largest eigenvalue, in absolute value, of a flat direction.

- tol:

  The largest decrement of a direction that leaves the verdict.

## Value

A list with `dirs`, the flat directions as vectors over every
coordinate, and `settled`, those whose own decrement is at most `tol`.

## Details

A coordinate the curvature does not resolve shows on the diagonal, which
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
reads one coordinate at a time. A direction along which several
coordinates move together can be flat while every coordinate alone is
curved, and the equilibrated matrix is what reads it: its diagonal is
one, so its smallest eigenvalue measures how far the coordinates are
from a linear dependence in the curvature, whatever their scales. A
direction is returned on the original coordinates, scaled so that its
largest entry is one in absolute value, with entries below 1e-3 set to
zero. Its own decrement is \\(q'\tilde g)^2/(2\|\lambda\|)\\, with \\q\\
the equilibrated eigenvector and \\\tilde g = g/s\\; the absolute value
is taken because the sign of an eigenvalue this small is the sign of
rounding. No direction is returned where an eigenvalue is below `-flat`:
the point is then not a maximum, and a flat direction there is not a
ridge of one.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`decrement_off()`](https://statmodels7.github.io/statmodels7/reference/decrement_off.md)
