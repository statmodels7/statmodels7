# Whether a Kinked Block's Curvature Is Read Centered

Returns, for an equation, whether its parametric intercept is profiled
out of the coordinate descent, so that a column's curvature is its
centered weighted sum of squares, as
[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md)
does for a descent.

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
derives this from the positions that a hold reaches through the
objective, and this function derives it from the `held_coef` of the
specification directly, because that slot is available in
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md).

## See also

[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md),
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
