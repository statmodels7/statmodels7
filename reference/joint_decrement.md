# The Newton Decrement of a Criterion at a Point

The Newton decrement \\\tfrac{1}{2}g^\top A^{-1} g\\, in the criterion's
own units.

## Usage

``` r
joint_decrement(g, A)
```

## Arguments

- g:

  The outer criterion's gradient, over the coordinates under test.

- A:

  The negated symmetric curvature over the same coordinates.

## Value

A single number, or `NA_real_`.

## Details

Under a quadratic model of the criterion at the reported point, taking
the step \\A^{-1}g\\ raises it by exactly this amount. The decrement is
therefore the part of the criterion that remains to be gained at the
point, in the units in which criteria are compared, and it approximates
the gap between the criterion at the point and at the optimum. A
gradient does not have that scale: its size depends on the number of
penalized coefficients and on the shape of the model, so that no single
threshold on it serves every model.

`NA` where \\A\\ is not positive definite, the point then not being a
maximum in those directions and no rise following from it.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`outer_curvature()`](https://statmodels7.github.io/statmodels7/reference/outer_curvature.md),
[`coord_decrement()`](https://statmodels7.github.io/statmodels7/reference/coord_decrement.md)
