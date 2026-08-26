# Is a Design Sparse, and the Zero Matrix to Accumulate It Into

`design_sparse()` reports whether any equation's block is sparse, and
`zero_information()` gives the square zero matrix of the right kind to
accumulate the information into.

## Usage

``` r
design_sparse(design)

as_dense(A)

as_sparse(A)

zero_information(design, total)
```

## Arguments

- design:

  The design, a list with one entry per distribution parameter, each
  carrying its block as `X`.

- total:

  The number of stacked coefficients across the equations, which is the
  side of the accumulator.

## Value

`design_sparse()` gives a single logical. `zero_information()` gives a
`total x total` matrix of zeros, a `dgCMatrix` when the design is sparse
and a base matrix otherwise.

## Details

The information is assembled one `crossprod` per parameter pair and each
product is written into a square accumulator. With a sparse design each
product is sparse, and writing a sparse block into a dense accumulator
signals that the number of items to replace is not a multiple of the
replacement length, from inside the assembly and naming nothing a caller
wrote. The accumulator therefore follows the design.

## See also

[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md),
which accumulates into it,
[`zap_nonfinite()`](https://statmodels7.github.io/statmodels7/reference/zap_nonfinite.md)
for the other place the two storages meet.
