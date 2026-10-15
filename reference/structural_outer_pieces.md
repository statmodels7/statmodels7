# The Penalty's Pieces on the Joint Vector

[`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)
for the vector that consists of the coefficients and the own parameters
of a filter. The determinant of a marginal criterion spans that vector
where a penalty covers those parameters.

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

The arithmetic is that of
[`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md);
the difference is where the coordinates of each unit are placed. An
ordinary unit is matched into the kept coefficients, a penalty over the
own parameters of a structural term is placed among the free parameters
of the filter, and a mixed covariance class is read where each of its
coordinates lives and put back in the interleaved order of the class.
These three readings are those of
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md),
and they are reused here so that the gradient and the Hessian place each
coordinate at the same position.

## See also

[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md),
[`outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/outer_pieces.md)
