# The Movement of the Predictors and the Information Along Directions

`jfp_eta_dir()` returns the movement of each equation's predictor along
a direction, the Jacobian times the direction; `jfp_dI()` a derivative
of the family's expected information, with some predictor indices fixed
and the others contracted against per-equation movements.

## Usage

``` r
jfp_eta_dir(st, v)

jfp_dI(st, q, r, fixed, dirs)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- v:

  A direction over the joint vector.

- q, r:

  The pair of the information.

- fixed:

  The equations of the derivatives that are not contracted.

- dirs:

  A list of movements, each a list of numeric vectors, one per
  distribution parameter.

## Value

`jfp_eta_dir()` a list of numeric vectors, one per distribution
parameter; `jfp_dI()` a numeric vector over the observations.
