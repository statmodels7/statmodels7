# A Basis of the Directions That a Set of Linear Constraints Leaves Free

The columns of a `p x r` matrix \\N\\ with \\B N = 0\\ and \\N^\top N =
I\\, spanning the null space of `B`.

## Usage

``` r
null_basis(B, p = ncol(B))
```

## Arguments

- B:

  A numeric matrix with `p` columns.

- p:

  The number of columns.

## Value

The matrix \\N\\, with zero columns where `B` has full column rank.

## Details

Where every row of `B` is a multiple of a unit vector, the constraints
hold coordinates, and \\N\\ is the identity's columns of the other
coordinates, so the reduced problem is the problem without those
coordinates exactly. Otherwise \\N\\ is read off the complete QR
decomposition of \\B^\top\\.
