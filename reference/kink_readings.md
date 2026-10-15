# The Inner Readings on the Directions That the Kinks Leave Free

[`certificate_readings()`](https://statmodels7.github.io/statmodels7/reference/certificate_readings.md)
over the null space of the linear constraints that a break-point sitting
on an observation imposes
([`kink_constraints()`](https://statmodels7.github.io/statmodels7/reference/kink_constraints.md)),
restricted to the coordinates in `keep`.

## Usage

``` r
kink_readings(g, A, keep = seq_along(g), kc = NULL)
```

## Arguments

- g:

  The score.

- A:

  The penalized information.

- keep:

  The coordinates read.

- kc:

  The constraint matrix, with as many columns as the coefficients
  (padded with zeros where `g` is longer), or `NULL`.

## Value

The two readings, as
[`certificate_readings()`](https://statmodels7.github.io/statmodels7/reference/certificate_readings.md)
returns them.

## Details

With no constraint this is
[`certificate_readings()`](https://statmodels7.github.io/statmodels7/reference/certificate_readings.md)
itself. A constraint that holds one coordinate removes it from the
reading. A constraint over a developed position, a row \\z_j^\top\\ of
its design, leaves a hyperplane, and the readings are those of the
reduced system \\(N^\top g, N^\top A N)\\ with \\N\\ an orthonormal
basis of the null space
([`null_basis()`](https://statmodels7.github.io/statmodels7/reference/null_basis.md)).
