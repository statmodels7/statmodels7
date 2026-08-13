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

  The design.

- total:

  The number of stacked coefficients.

## Value

A logical, or a square matrix of zeros.

## Details

The information is assembled one `crossprod` per parameter pair and
placed into a square accumulator. With a sparse design each product is
sparse, and placing it into a dense accumulator signals that the number
of items to replace is not a multiple of the replacement length – from
inside the assembly, naming nothing a caller wrote. The accumulator
follows the design instead.
