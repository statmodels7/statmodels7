# Store a Small Design as Base Matrices

Converts the sparse blocks of a design to base matrices when the
equations together carry fewer than `min_dim` coefficients and no block
is sparse because the caller requested it.

## Usage

``` r
densify_small(design, min_dim = 100L)
```

## Arguments

- design:

  The design, a list with one entry per distribution parameter, each
  carrying its block as `X`.

- min_dim:

  The number of coefficients below which the design is stored dense.

## Value

The design, with every `X` a base matrix where the total is below
`min_dim`, and unchanged otherwise.

## Details

A random effect builds its block sparse whatever its size. Below about a
hundred coefficients the sparse route costs more than it saves, its
fixed cost being the coercions and the S4 dispatch around each product,
which do not shrink with the matrix;
[`worth_sparse()`](https://statmodels7.github.io/statmodels7/reference/worth_sparse.md)
records the same crossover for the factorization. The conversion is
applied only where no block moves with the coefficients and no
structural term is present, since those designs are rebuilt during the
fit and a rebuilt block would come back sparse. A block that the caller
requested sparse (`sparse = TRUE` on a term,
`linpar_options(sparse = TRUE)`, a sparse matrix as input) is never
converted, and then neither is the rest of the design: the design keeps
the storage it was given.
