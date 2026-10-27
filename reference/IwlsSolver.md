# S7 Classes for the Solver of a Scoring Step

`IwlsSolver` is the abstract parent of the objects that select how the
linear system of one step of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
is solved. Its three subclasses are `SolveQR`, `SolveSVD` and
`SolveChol`. They are built by
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md),
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md),
which are the functions to call; the classes are exported so that
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html)
can test an object.
[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
compares the three.

`SolveQR` has the properties `tol` and `floor`, `SolveSVD` has `tol` and
`floor`, and `SolveChol` has `floor`. The raw constructors do not
validate them; the functions
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md),
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
do.

## Usage

``` r
IwlsSolver()

SolveQR(tol = integer(0), floor = integer(0))

SolveSVD(tol = NULL, floor = integer(0))

SolveChol(floor = integer(0))
```

## Arguments

- tol:

  For `SolveQR`, the tolerance of the pivoted QR, a single positive
  number. For `SolveSVD`, `NULL` or the relative cut of the singular
  values, a single positive number.

- floor:

  The relative floor of the eigenvalues of a penalized matrix that is
  not positive definite, a single positive number.

## Value

An S7 object of the class named, which inherits from `IwlsSolver`.

## See also

[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
for the comparison of the three;
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md),
[`solve_svd()`](https://statmodels7.github.io/statmodels7/reference/solve_svd.md)
and
[`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)
for the constructors.

## Examples

``` r
S7::S7_inherits(solve_svd(), IwlsSolver)
#> [1] TRUE
S7::S7_inherits(solve_svd(), SolveSVD)
#> [1] TRUE
solve_svd()@tol
#> NULL
```
