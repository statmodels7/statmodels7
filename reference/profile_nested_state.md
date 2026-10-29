# The Nested Problem at a Point

Prepares what the expansion of the quantities that the outer criterion
estimates needs at a point of the profile. The held quantities are the
coefficients and the hyperparameters of `which`; the unknown ones are
the coefficients that the criterion estimates and, with
`hyper = "refit"`, the other hyperparameters. Both sets are coordinates
of the criterion, and the function holds the exact Hessian of the
criterion over them, with the functions that solve the inner problem and
read the exact gradient of the criterion at given values.

## Usage

``` r
profile_nested_state(ctx, pt)
```

## Arguments

- ctx:

  The context of the call.

- pt:

  A point, as
  [`profile_point()`](https://statmodels7.github.io/statmodels7/reference/profile_point.md)
  returns it, with the held values in `val`.

## Value

A list, or `NULL` where the Hessian of the criterion is not available.

## Details

With the held quantities \\h\\ and the unknown ones \\u\\, the criterion
\\V(u, h)\\ is evaluated at the penalized mode of the other
coefficients, and none of the coefficient coordinates enters its
determinant. The result of
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
over all of them holds the blocks \\V\_{uu}\\ and \\V\_{uh}\\, and the
derivative of the unknown quantities is \\du/dh = -V\_{uu}^{-1}
V\_{uh}\\. The hyperparameters are on the scale of their links, the
scale of the grid.
