# The Variance of the Quantities a Fit Reports

The delta method over the joint variance, with the Jacobian
[`readable_joint`](https://statmodels7.github.io/statmodels7/reference/readable_joint.md)
supplies.

## Usage

``` r
readable_vcov(spec, design, fit, V)
```

## Arguments

- spec:

  The fitted specification.

- design:

  The design.

- fit:

  The fit.

- V:

  The joint variance over the coordinates.

## Value

The variance over the quantities.

## Details

A quantity that reads a coordinate whose variance is missing – one a
kinked penalty set to zero, or a parameter held by an intercept – has no
variance either, and its row and column are left missing rather than
computed from the coordinates that do have one.
