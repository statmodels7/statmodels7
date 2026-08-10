# Cholesky Factors of Small Blocks, Vectorized Over Observations

Returns the lower triangular \\L_i\\ with \\L_iL_i' = \Omega_i\\ for
every observation, as a \\n \times K \times K\\ array.

## Usage

``` r
chol_blocks(Om)
```

## Arguments

- Om:

  An \\n \times K \times K\\ array of symmetric blocks.

## Value

An array of the same shape, lower triangular in its last two indices, or
`NULL` when some block is not positive definite.

## Details

The standard recursion is written over the \\K\\ indices, which are few,
and evaluated over all observations at once, which are many.

A block that is not positive definite has a non-positive pivot, and the
function returns `NULL` at that point rather than taking its square
root: the observed curvature far from the optimum is routinely
indefinite, so this is an ordinary outcome the caller answers by falling
back to the assembled route, and a warning about a `NaN` would report it
as a defect.
