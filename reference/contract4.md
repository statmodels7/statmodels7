# The Fourth Derivative of the Objective Contracted Twice

\\U\[v, u\] = (\partial^2 K/\partial\beta^2)\cdot(v, u)\\, a matrix over
the stacked coefficients.

## Usage

``` r
contract4(spec, design, d4, params, npar, offs, total, tv, tu)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- d4:

  The fourth derivatives in the link scale.

- params, npar, offs, total:

  The block bookkeeping.

- tv, tu:

  The predictors of the two directions.

## Value

A square matrix.
