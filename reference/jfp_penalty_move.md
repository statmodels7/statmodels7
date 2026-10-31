# How the Jeffreys Curvature Moves With the Coefficients

Returns the derivative of the Jeffreys penalty's Hessian along a
direction \\v\\, or its second derivative along \\v\\ and \\w\\, with
the sign of a penalty: \\-\sum_b a_b D_3(v)\\ and \\-\sum_b a_b D_4(v,
w)\\ of
[`jfp_deriv3_dir()`](https://statmodels7.github.io/statmodels7/reference/jfp_deriv3_dir.md).
With `M` the trace against it is returned instead of the matrix.

## Usage

``` r
jfp_penalty_move(spec, design, coef, v, w = NULL, M = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  A named list of coefficient vectors.

- v, w:

  Directions over the coefficients; `w = NULL` for the first derivative.

- M:

  A symmetric matrix over the coefficients, or `NULL`.

## Value

A base matrix, or a single number when `M` is given; `NULL` (or `0` when
`M` is given) where the model carries no Jeffreys term.

## See also

[`statmod_penalty_dbeta()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_dbeta.md),
[`statmod_penalty_second()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_second.md)
