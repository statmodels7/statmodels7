# Print a Solver

Writes one line naming the factorization and its tolerances.

## Usage

``` r
# S3 method for class 'IwlsSolver'
print(x, ...)
```

## Arguments

- x:

  An object inheriting from
  [IwlsSolver](https://statmodels7.github.io/statmodels7/reference/IwlsSolver.md).

- ...:

  Unused.

## Value

`x`, invisibly. The function is called for the output that it writes.

## See also

[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)

## Examples

``` r
solve_qr()
#> QR solver (tol 1e-07, floor 1e-08)
```
