# The Penalty's Pieces on the Joint Vector

[`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)
over the coefficients AND a filter's own parameters, which is the vector
a marginal criterion's determinant spans where a penalty covers those
parameters.

## Usage

``` r
structural_outer_pieces(spec, design, coef, hyper, idx, jd)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- idx:

  The outer index.

- jd:

  The joint rows.

## Value

A list with `S`, `c`, `S2`, `c2`, `rho2` and `pair`, in the shape
[`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)
returns them.

## Details

The arithmetic is
[`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)'s
and the difference is where each unit's coordinates live: an ordinary
unit is matched into the kept coefficients, a penalty over a structural
term's own parameters is placed among the filter's free ones, and a
MIXED covariance class is read where each of its coordinates lives and
put back in the class's own interleaved order. Those three readings are
the gradient's own, in
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md),
and are composed here so the two cannot disagree about a position.

## See also

[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md),
[`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)
