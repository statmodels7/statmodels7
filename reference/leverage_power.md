# Apply a Function of a Symmetric Block Through Its Eigenvalues

Computes \\(I - H)^{p}\\ for a symmetric \\H\\, with the factor zero on
the directions where an eigenvalue of \\H\\ reaches one.

## Usage

``` r
leverage_power(H, power)
```

## Arguments

- H:

  A symmetric matrix.

- power:

  The exponent, negative.

## Value

A list with the matrix `F` and `capped`, the number of eigenvalues
treated as one.
