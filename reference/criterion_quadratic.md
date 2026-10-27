# The Quadratic Expansion of the Outer Criterion

Reads the exact outer gradient and the exact outer Hessian of a marginal
criterion at the fit and returns the entries that belong to one
hyperparameter.

## Usage

``` r
criterion_quadratic(fit, cd)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- cd:

  The list returned by
  [`criterion_setup()`](https://statmodels7.github.io/statmodels7/reference/criterion_setup.md).

## Value

A list with the gradient `g` and the curvature `A` (minus the Hessian),
both on the free coordinate of the hyperparameter.
