# The Hyperparameters an Outer Method Estimates

One row per hyperparameter of a penalized term that is fitted in the
joint system, with the link that carries it onto the whole line.

## Usage

``` r
outer_hyper_index(spec, blocks)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- blocks:

  The block split, as
  [`statmod_blocks`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns it.

## Value

A data frame with `parameter`, `term` and `name`, and a list column-free
`link` carried as an attribute list.

## Details

A term whose penalty has a kink is left out. Its coefficients are
estimated by a proximal method with everything else held fixed, and the
criterion is a Laplace approximation, which asks for a second derivative
that does not exist at a kink.
