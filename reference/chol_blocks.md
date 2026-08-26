# Cholesky Factors of Small Blocks, Vectorized Over Observations

Returns the lower triangular \\L_i\\ with \\L_iL_i' = \Omega_i\\ for
every observation, as a \\n \times K \times K\\ array.

## Usage

``` r
chol_blocks(Om)
```

## Arguments

- Om:

  An \\n \times K \times K\\ array of symmetric blocks, as
  [`info_blocks()`](https://statmodels7.github.io/statmodels7/reference/info_blocks.md)
  returns it. Only the lower triangle of each block is read.

## Value

An \\n \times K \times K\\ numeric array, lower triangular in its last
two indices, with `Om[i, , ] == L[i, , ] %*% t(L[i, , ])` for every `i`.
`NULL` as soon as any block fails to be positive definite, so a single
bad observation declines for the whole sample.

## Details

\\K\\ is the number of distribution parameters, so it is 1, 2 or 3 for
most families and never large; \\n\\ is the number of observations and
is. The standard Cholesky recursion is therefore written out over the
\\K\\ indices and evaluated over all \\n\\ observations at once, one
vectorized pass per entry of the factor, so no loop runs over
observations.

A block that is not positive definite produces a non-positive pivot, and
the function returns `NULL` at that point instead of taking its square
root. That is an expected outcome: the observed curvature far from the
optimum is routinely indefinite, and the caller answers by falling back
to the assembled route. A warning about a `NaN` would report an ordinary
branch as a defect.

## See also

[`info_blocks()`](https://statmodels7.github.io/statmodels7/reference/info_blocks.md)
for the input,
[`sqrt_design()`](https://statmodels7.github.io/statmodels7/reference/sqrt_design.md)
for what the factors are used to build.
