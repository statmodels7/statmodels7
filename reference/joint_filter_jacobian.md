# The Derivative Rows of Every Equation Under a Joint Filter

Returns the derivative of each equation's predictor in the coefficients
of every equation and the joint filter's free parameters, at the fitted
point, together with the first derivative of the driving vector in every
predictor.

## Usage

``` r
joint_filter_jacobian(spec, design, coef)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

## Value

`NULL` where the model carries no joint filter, otherwise a list with
`V` (one matrix per distribution parameter), `U1` (the derivative of
each driving score in every predictor), `G` (the positions of the
filtered parameters), `nb`, `key` (the names of the columns in the
variance matrix), `free` (the filter's free parameters), `f` and `ev`.

## Details

The rows of a filtered equation are the forward Jacobian of the joint
recursion from
[`joint_filter_curvature()`](https://statmodels7.github.io/statmodels7/reference/joint_filter_curvature.md);
those of the other equations are their static designs. The columns of a
held level are left out, as the variance matrix leaves them out.

## See also

[`structural_se_columns()`](https://statmodels7.github.io/statmodels7/reference/structural_se_columns.md),
[`joint_continued_deriv_inputs()`](https://statmodels7.github.io/statmodels7/reference/joint_continued_deriv_inputs.md)
