# The Pieces From Which the Terms Are Built

`jfp_eta_dir()` returns the movement of each equation's predictor along
a direction; `jfp_dI()` a derivative of the family's expected
information, contracted along directions; `jfp_white()`, `jfp_rows_B()`
and `jfp_pair_weights()` the three shapes in which a derivative of the
block's matrix enters (see the comment at the top of `R/jfp_outer.R`).

## Usage

``` r
jfp_eta_dir(st, v)

jfp_dI(st, q, r, fixed, dirs)

jfp_white(st, bs, dirs)

jfp_rows_B(st, bs, dirs)

jfp_pair_weights(st, bs, Y, dirs)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- v:

  A direction over the joint vector.

- q, r:

  The pair of the information.

- fixed:

  The equations of the derivatives that are not contracted.

- dirs:

  A list of directions, each as `jfp_eta_dir()` returns it.

- bs:

  A block's entry of the state.

- Y:

  A symmetric \\\tilde p\times\tilde p\\ matrix.

## Value

`jfp_eta_dir()` a list of numeric vectors, one per distribution
parameter; `jfp_dI()` a numeric vector over the observations;
`jfp_white()` a \\\tilde p\times\tilde p\\ matrix; `jfp_rows_B()` a
matrix over the joint vector and the \\\tilde p^2\\ entries;
`jfp_pair_weights()` a list keyed as
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
reads its `H`.
