# The Driving Vector's Derivatives in Every Predictor, at the Fit

Returns the data that
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
and its siblings read for a joint filter: the driving vector at the
fitted predictors and its derivatives in the predictors of every
equation, to the order that `order` gives.

## Usage

``` r
joint_drive_data(
  spec,
  tm,
  gl,
  H,
  D3 = NULL,
  D4 = NULL,
  D5 = NULL,
  order = 2L,
  theta = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- tm:

  The joint filter.

- gl, H, D3, D4, D5:

  The family's link-scale derivatives at the fitted predictors; those
  above the order that `order` gives may be `NULL`.

- order:

  The highest order of the predictor's derivative to be propagated: 1
  for the reverse pass, 2 for the curvature, 3 for the third derivative,
  4 for the fourth.

- theta:

  The per-observation parameters at the fitted predictors, from which
  the expected information is read where `scaling` is not zero.

## Value

A list with `filtered`, `u`, `U1`, `U2` and, to the order given, `U3`
and `U4`.

## Details

With `scaling` zero the driving vector is the link-scale gradient of the
filtered parameters, and its derivatives are the family's own: the entry
\\U\_{l;a_1\dots a_r}\\ is the derivative of order \\r + 1\\ with
indices \\l, a_1, \dots, a_r\\. Otherwise
[`joint_drive_scaled()`](https://statmodels7.github.io/statmodels7/reference/joint_drive_scaled.md)
composes them with the derivatives of the scaling matrix.
