# The Pieces of the Block Variance Over the Estimated Coefficients

The quantities
[`marginal_vcov()`](https://statmodels7.github.io/statmodels7/reference/marginal_vcov.md)
and
[`marginal_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/marginal_edf_correction.md)
share: the inverse of the penalized information over the integrated
coefficients, the movement of those coefficients with the estimated ones
and with the hyperparameters, and the Hessian of the criterion over
\\(\eta, \gamma)\\.

## Usage

``` r
marginal_blocks(mc, coef, hyper, method, A, keep, full = TRUE)
```

## Arguments

- mc:

  The list
  [`fit_marginal_context()`](https://statmodels7.github.io/statmodels7/reference/fit_marginal_context.md)
  returns.

- coef, hyper:

  The fitted coefficients and hyperparameters.

- method:

  The
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  the fit ran.

- A:

  The penalized information over the kept coordinates.

- keep:

  A logical vector over every coefficient, `TRUE` where it is kept.

- full:

  Whether to build `Be`, the movement with the hyperparameters.

## Value

A list with `Pu`, `ui`, `p`, `Bg`, `Be`, `Hv`, `nh`, `ng` and `idx`, or
`NULL` where a piece cannot be computed.

## Details

Write \\u\\ for the coefficients the criterion integrates and \\\gamma\\
for the ones it estimates. At the mode, \\\partial u/\partial\gamma =
-A\_{uu}^{-1}A\_{u\gamma}\\, and a hyperparameter moves \\u\\ by
\\-A\_{uu}^{-1}\\\partial^2\rho/\partial u\\\partial\eta\\. The columns
of `Bg` and `Be` are these movements written over the kept coordinates,
with a unit entry at \\\gamma\\'s own position.
