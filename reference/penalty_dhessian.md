# The Derivative of a Penalty's Hessian in Its Hyperparameters

One matrix per hyperparameter, for a penalty whose Hessian is linear in
them.

## Usage

``` r
penalty_dhessian(pen, beta, theta)
```

## Arguments

- pen:

  A penalties7 penalty.

- beta:

  The term's coefficients.

- theta:

  The term's hyperparameters.

## Value

A named list of matrices.

## Details

Under that linearity \\S(\theta + h e_m) = S(\theta) + h\\\partial
S/\partial\theta_m\\ holds for every \\h\\, so the difference below is
exact arithmetic and not an approximation: there is no truncation error
to shrink and no step to choose. It is written as a difference rather
than by evaluating at a unit vector because a hyperparameter of zero is
outside the bounds a penalty validates against, while \\2\theta_m\\ is
not.

The linearity is verified by
[`penalty_hessian_linear`](https://statmodels7.github.io/statmodels7/reference/penalty_hessian_linear.md)
before this is called.
