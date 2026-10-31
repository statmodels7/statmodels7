# A Solve With a Matrix Plus a Product of Low Rank

`lowrank_solve()` solves \\(A + VV^\top + \lambda I)\delta = u\\ by
Woodbury's identity, factoring \\A + \lambda I\\ in its own storage and
a small matrix of the size of \\V\\'s columns. `lowrank_dense()` returns
\\A + VV^\top\\ as a base matrix.

## Usage

``` r
lowrank_solve(A, V, u, damp = 0)

lowrank_dense(A, V)
```

## Arguments

- A:

  A symmetric matrix, dense or a Matrix object.

- V:

  A matrix with as many rows as `A`, or `NULL`.

- u:

  The right-hand side.

- damp:

  The Levenberg damping added to the diagonal.

## Value

`lowrank_solve()`: a numeric vector, or `NULL`. `lowrank_dense()`: a
base matrix.

## Details

With \\A\\ positive definite, \$\$(A + VV^\top)^{-1}u = A^{-1}u -
A^{-1}V(I + V^\top A^{-1}V)^{-1} V^\top A^{-1}u,\$\$ and the matrix in
the middle is positive definite. The Jeffreys curvature enters the
assembled route in this form
([`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)),
so a sparse information keeps its sparse factorization. Where \\A\\ has
no Cholesky factor the function returns `NULL`, and
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
assembles the sum and repairs it as it repairs any indefinite matrix.

## See also

[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
