# The Newton Decrement of a Joint Fit, With the Verdict on Its Curvature

Reads, at a point, how far a joint objective over the coefficients and a
structural term's parameters lies above its maximum, and whether the
curvature there admits such a reading.

## Usage

``` r
joint_reading(g, K, keep, own = integer(0))
```

## Arguments

- g:

  The gradient of the joint objective.

- K:

  Its curvature, a matrix as long as `g`.

- keep:

  The coordinates to read.

- own:

  The coordinates, among those of `g`, that are the structural term's
  own parameters.

## Value

A list with `decrement` (`NA` where it cannot be computed, when `ascent`
is `TRUE`, and when `partial` is `TRUE` and the lower bound is below
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)),
`ascent`, `ratio`, the smallest eigenvalue's magnitude in units of
\\\varepsilon\lambda\_{\max}\\ (`NA` when there is nothing to read), and
`partial`.

## Details

The reading is the Newton decrement \\\tfrac12 g^\top K^{-1} g\\ over
the coordinates in `keep`, after three steps taken in this order.

A coordinate whose diagonal of \\K\\ is not positive is flat at an edge,
with the sign set by rounding, and is left out. If no coordinate has a
positive diagonal the curvature is not flat but negative, nothing is
left to read, and `ascent` is `TRUE`.

The remaining curvature is read for a direction of ascent. If its
smallest eigenvalue is below \\-c\\\varepsilon\\\lambda\_{\max}\\, with
\\c\\ from
[`curvature_ascent_limit()`](https://statmodels7.github.io/statmodels7/reference/curvature_ascent_limit.md),
the point is not a maximum, no decrement follows from an indefinite
matrix, and `ascent` is `TRUE`.
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
is not consulted in that case, because on an indefinite matrix it can
name the term's own parameters, and the decrement over the coordinates
that remain is then smaller than the distance to the maximum.

Otherwise the matrix is positive semidefinite up to rounding, the
coordinates that
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
names are left out as aliased, and the decrement is computed over the
others. A negative decrement is not possible for a positive semidefinite
matrix, so it is read as a direction of ascent.

A column of the design that duplicates another is aliased, and its
coordinate carries no information because the objective is flat along
it. The same holds for no parameter of the term: the curvature of a
filter that has exploded is numerically of rank one, and
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
then names the term's own parameters, which the objective is not flat
in. Where it names one of the coordinates in `own`, the decrement over
the others is a lower bound of the decrement over all of them, so it can
show that the point is not a maximum and cannot show that it is one: a
value above
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
is returned and a value below it is replaced by `NA`, with `partial` set
to `TRUE`.

## See also

[`curvature_ascent_limit()`](https://statmodels7.github.io/statmodels7/reference/curvature_ascent_limit.md),
[`alternation_readings()`](https://statmodels7.github.io/statmodels7/reference/alternation_readings.md),
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`filter_decrement()`](https://statmodels7.github.io/statmodels7/reference/filter_decrement.md)
