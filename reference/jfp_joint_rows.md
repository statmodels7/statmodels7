# The Derivative of Every Predictor Beside a Filter

Returns, for a model carrying a filter, the derivative of each
equation's predictor in the joint vector of the coefficients and the
filter's free parameters, with what the second derivative of a filtered
predictor needs.

## Usage

``` r
jfp_joint_rows(spec, design, coef, ev)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design as it arrived, carrying the structural state.

- coef:

  The coefficients.

- ev:

  What
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns at `coef`.

## Value

A list with `V` (one dense matrix per distribution parameter over the
kept coordinates), `keep`, `driven` (the equations whose predictor the
filter writes), and `curv`, a function of the per-observation weights (a
named list keyed by parameter) returning \\\sum_t\sum_c
g\_{c,t}\\\partial^2\eta\_{c,t}/\partial u\\\partial u\\ over the kept
coordinates.

## Details

The rows are those of
[`joint_design_rows()`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md),
with the equation that the filter drives replaced by the forward
Jacobian of the recursion
([`filter_curvature()`](https://statmodels7.github.io/statmodels7/reference/filter_curvature.md),
or
[`joint_filter_curvature()`](https://statmodels7.github.io/statmodels7/reference/joint_filter_curvature.md)
for a joint filter), and a level that an intercept holds dropped as
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
drops it.
