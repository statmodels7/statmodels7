# The Nested Problem at a Point

Prepares what the expansion of the coordinates that the outer criterion
estimates needs at a point of the profile: the held coefficients and
those coordinates as one set of criterion coordinates, the exact Hessian
of the criterion over that set, and the functions that solve the inner
problem and read the exact gradient of the criterion at given values.

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

With the held coefficients \\\theta\\ and the coordinates \\\gamma\\
that the criterion estimates, the criterion \\V(\gamma, \theta)\\ is
evaluated at the penalized mode of the other coefficients, and none of
the two sets enters its determinant. Both sets are therefore criterion
coordinates of
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
whose result over \\(\gamma, \theta)\\ holds the blocks
\\V\_{\gamma\gamma}\\ and \\V\_{\gamma\theta}\\. The derivative of the
estimated coordinates is \\d\gamma/d\theta = -V\_{\gamma\gamma}^{-1}
V\_{\gamma\theta}\\.
