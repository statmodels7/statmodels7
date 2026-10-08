# The Forward Jacobian of a Filter in Every Coordinate

The derivative of the predictor that a score-driven term produces, in
the coefficients of every equation followed by the term's free
parameters on their unconstrained scale, at each observation, together
with the static rows of every equation in the same columns.

## Usage

``` r
filter_joint_jacobian(spec, design, coef)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- coef:

  The coefficients.

## Value

`NULL` where the model carries no filter, otherwise a list with `J` (the
filter equation's rows), `V` (every equation's rows, the filter's being
`J`), `H` (the family's second derivatives on the link scale), `ap` (the
index of the filter's parameter), `nb` (the number of coefficient
columns), `key` (the names of the columns in the variance matrix) and
`free` (the term's free parameters), all columns restricted to the free
ones.

## Details

The rows are the ones the joint information is assembled from
([`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)):
each equation's design placed in its own columns, and for the filter's
equation the Jacobian that
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
propagates beside the state. A model carries at most one structural
term, so every other equation is static and its design is its
derivative.

## See also

[`structural_se_columns()`](https://statmodels7.github.io/statmodels7/reference/structural_se_columns.md),
[`continued_deriv_inputs()`](https://statmodels7.github.io/statmodels7/reference/continued_deriv_inputs.md)
