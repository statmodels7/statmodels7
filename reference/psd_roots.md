# Symmetric Square Roots of a Stack of Small Information Matrices

For each slice of an array of symmetric matrices, computes the symmetric
square root and its pseudo-inverse from the eigen decomposition.

## Usage

``` r
psd_roots(W)
```

## Arguments

- W:

  An array of dimension \\n \times d \times d\\.

## Value

A list with `root` and `root_pinv`, arrays of the same dimension, and
`clipped`, the number of slices with a negative eigenvalue.

## Details

Eigenvalues below \\10^{-12}\\ times the largest eigenvalue of the slice
are set to zero and have a zero pseudo-inverse. This covers the negative
eigenvalues of an observed information and the zero information of an
observation with weight zero.
