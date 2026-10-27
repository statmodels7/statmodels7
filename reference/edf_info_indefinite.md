# Whether an Information Matrix Is Indefinite

Returns `TRUE` when the symmetric matrix `H` has an eigenvalue below
\\-p\\\epsilon\\\lambda\_{\max}\\, with \\p\\ the dimension and
\\\epsilon\\ the machine epsilon. A matrix with a Cholesky factor is
positive definite and is not decomposed.

## Usage

``` r
edf_info_indefinite(H)
```

## Arguments

- H:

  A square numeric matrix.

## Value

A single logical; `FALSE` where `H` has a non-finite entry.
