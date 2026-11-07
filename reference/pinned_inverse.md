# The Inverse of a Penalized Information with Pinned Coordinates

Inverts a penalized information after pinning the given coordinates, and
those found deficient, to the identity, and returns the inverse with the
pinned rows and columns set to zero.

## Usage

``` r
pinned_inverse(K, hold = integer(0))
```

## Arguments

- K:

  The penalized information, a square matrix.

- hold:

  The coordinates to pin.

## Value

A matrix of the dimension of `K`, or `NULL` if no factorization exists
after the deficient coordinates are pinned.

## Details

The coordinates are pinned by
[`pin_boundary()`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md),
as
[`ctx_penalized()`](https://statmodels7.github.io/statmodels7/reference/ctx_penalized.md)
pins the matrix that governs the movement of the mode. If the Cholesky
factorization fails, the coordinates that
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
reports are pinned as well and the factorization is repeated. The zero
rows and columns state that the mode does not move in a pinned
coordinate, as in
[`mode_solver()`](https://statmodels7.github.io/statmodels7/reference/mode_solver.md).

## See also

[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md),
the caller.
