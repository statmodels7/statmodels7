# Carry a Hyperparameter Derivative onto the Free Scale

The diagonal chain rule each hyperparameter's own link induces, at first
order for a gradient and at second order for a Hessian.

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

A vector, or a matrix.
