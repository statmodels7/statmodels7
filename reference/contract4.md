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

## Details

The blocks follow
[`contract3()`](https://statmodels7.github.io/statmodels7/reference/contract3.md):
each is minus a weighted crossproduct of the two designs, with the
weight formed from the fourth derivatives and the predictors of the two
directions.
