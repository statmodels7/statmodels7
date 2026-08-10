# The Per-Hyperparameter Pieces of the Outer Derivatives

The penalty's contributions to the criterion's gradient and Hessian,
placed in the stacked coefficient space.

## Usage

``` r
outer_pieces(spec, design, coef, hyper, idx, offs, total, order = 2L)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- idx:

  The outer index.

- offs, total:

  The block offsets and the total width.

- order:

  `1` for the first-order pieces alone, `2` for the second-order ones as
  well.

## Value

A list with `S` (one matrix per hyperparameter) and `c` (one vector),
and at order 2 also `S2` and `c2` (one per pair, keyed), `rho2` (the
hyperparameter Hessian) and `pair` (the key of each pair).

## Details

At order 1 the second-order generics are not called at all. That is what
lets a penalty supplying only `penalty_dhessian()` give an exact
gradient: asking it for a derivative the gradient does not use would
have rejected it for a quantity nobody wanted.
