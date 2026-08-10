# The Per-Hyperparameter Pieces of the Outer Derivatives

The penalty's contributions to the criterion's gradient and Hessian,
placed in the stacked coefficient space.

## Usage

``` r
outer_pieces(spec, design, coef, hyper, idx, offs, total)
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

## Value

A list with `S` (one matrix per hyperparameter), `c` (one vector), `S2`
and `c2` (one per pair, keyed), `rho2` (the hyperparameter Hessian) and
`pair` (the key of each pair).
