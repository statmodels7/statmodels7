# Run the Compiled Coordinate Descent on Either Storage

Sends the block to the dense kernel or to the sparse one, taking a
`dgCMatrix` apart into the slots the second reads.

## Usage

``` r
coord_call(X, z, w, b0, tab, screen, tol, covariance)
```

## Arguments

- X:

  The block, `n x p`, dense or `dgCMatrix`.

- z, w:

  The working response and weights, each of length `n`.

- b0:

  The starting coefficients, length `p`.

- tab:

  The piecewise linear proximal table, as
  [`penalties7::penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.html)
  returns it.

- screen:

  The **zero-based** positions the strong rule kept, for the C++
  indexing.

- tol:

  The stopping tolerance on the largest coefficient change of a sweep.

- covariance:

  `TRUE` to hold the gradient itself and cache Gram columns, `FALSE` to
  hold the running residual.
  [`coord_covariance()`](https://statmodels7.github.io/statmodels7/reference/coord_covariance.md)
  decides.

## Value

The kernel's list of three: `beta` (the fitted coefficients, length
`p`), `sweeps` (how many passes it took) and `grad` (the gradient at the
point reached, length `p`).

## Details

The two kernels are one algorithm instantiated twice over a column
accessor, and they agree bit for bit. The arithmetic licenses that:
skipping a structural zero omits an addition of zero, which is exact. It
is the one place in this toolkit where an identity assertion over
compiled floating point is correct.

The `dgCMatrix` is taken apart in R, so the compiled code needs no
dependency on the Matrix package's C API.

## See also

[`coord_block()`](https://statmodels7.github.io/statmodels7/reference/coord_block.md)
