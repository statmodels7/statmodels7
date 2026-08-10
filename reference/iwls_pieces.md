# The Pieces One Scoring Step Needs

Assembles, at the current coefficients, whichever of the square-root
design, the penalty's factor and the penalized information the requested
decomposition will actually use.

## Usage

``` r
iwls_pieces(spec, design, coef, hyper, method)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients, per parameter.

- hyper:

  The hyperparameters.

- method:

  An
  [`Iwls`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

## Value

A list with `R`, `C` and `A`.
