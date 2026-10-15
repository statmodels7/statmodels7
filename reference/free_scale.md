# Carry a Hyperparameter Derivative onto the Free Scale

The diagonal chain rule that each hyperparameter's own link induces, at
first order for a gradient and at second order for a Hessian.

## Usage

``` r
free_scale(g, hyper, idx)

free_scale2(M, g, hyper, idx)

link_slopes(hyper, idx)
```

## Arguments

- g:

  A gradient on the parameter scale.

- hyper:

  The hyperparameters.

- idx:

  The outer index.

- M:

  A Hessian on the parameter scale.

## Value

For `free_scale()` a vector, for `free_scale2()` a matrix, and for
`link_slopes()` a list with the components `h1` and `h2`, the first and
second derivatives of each hyperparameter's inverse link at its value.
