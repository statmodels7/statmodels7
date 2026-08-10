# Is a Penalty's Hessian Linear in the Hyperparameters?

Asks the penalty rather than assuming: its Hessian must double when the
hyperparameters double, and must not move when the coefficients do.

## Usage

``` r
penalty_hessian_linear(pen)
```

## Arguments

- pen:

  A penalties7 penalty.

## Value

A single logical.
