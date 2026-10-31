# The Criterion's Determinant With a Jeffreys Curvature

Computes \\\log\|K\|\\ and the dimension \\q\\ that
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
reads, for a model carrying a Jeffreys term, with \\K = K_0 + VV^\top\\
factored in its two parts.

## Usage

``` r
marginal_lowrank(spec, design, coef, hyper, H, basis = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef, hyper:

  The coefficients and the hyperparameters.

- H:

  The information.

- basis:

  The integrated basis, or `NULL`.

## Value

A list with `logdet` and `q`, or `NULL` where the model carries no
Jeffreys term or \\K_0\\ has no positive definite factor.

## Details

\\\log\|K\| = \log\|K_0\| + \log\|I + V^\top K_0^{-1}V\|\\, the matrix
determinant lemma. Under
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
both parts are projected on the integrated basis first. The coordinates
are pinned as
[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
pins them, and the pinned rows of \\V\\ are zero.

## See also

[`penalized_lowrank()`](https://statmodels7.github.io/statmodels7/reference/penalized_lowrank.md)
