# Whether a Kinked Block's Curvature Is Read Centered

Answers, for an equation, the question
[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md)
answers for a coordinate descent: whether its parametric intercept is
profiled out of the descent, so that a column's curvature is its
centered weighted sum of squares.

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
asks the same of the positions a hold reaches through the objective;
this asks it of the specification's `held_coef` directly, which is what
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
has in hand.

## See also

[`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md),
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
