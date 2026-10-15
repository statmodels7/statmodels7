# Whether a Kinked Block's Curvature Is Read Centered

Decides, for an equation, whether its parametric intercept is profiled
out of the coordinate descent, so that a column's curvature is its
centered weighted sum of squares, as
[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md)
decides for a descent.

## Usage

``` r
curv_centered(spec, design, p)
```

## Arguments

- spec, design:

  The specification and its design.

- p:

  The equation, a parameter name.

## Value

`TRUE` or `FALSE`.

## Details

The intercept is profiled out where the equation has one and it is not
held at a value.
[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md)
decides this from the positions that a hold reaches through the
objective; this function decides it from the `held_coef` of the
specification directly, which is what
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
has in hand.

## See also

[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md),
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
