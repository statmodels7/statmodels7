# Run the Compiled Coordinate Descent on Either Storage

Sends the block to the dense kernel or to the sparse one, taking a
`dgCMatrix` apart into the slots the second reads.

## Usage

``` r
coord_call(X, z, w, b0, tab, screen, tol, covariance)
```

## Arguments

- X:

  The block, dense or `dgCMatrix`.

- z, w:

  The working response and weights.

- b0:

  The starting coefficients.

- tab:

  The proximal table, from
  [`penalty_prox_spec`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.html).

- screen:

  The zero-based positions the strong rule kept.

- tol:

  The stopping tolerance on the coefficient change.

- covariance:

  Whether to hold the gradient rather than the residual.

## Value

The kernel's list: `beta`, `sweeps`, `grad`.

## Details

The two kernels are one algorithm instantiated twice over a column
accessor, so they agree BIT FOR BIT rather than to a tolerance: skipping
a structural zero omits an addition of zero, which is exact. The
`dgCMatrix` is decomposed here rather than in C++ so that the compiled
code needs no dependency on the Matrix package's C API.

## See also

[`coord_block`](https://statmodels7.github.io/statmodels7/reference/coord_block.md)
