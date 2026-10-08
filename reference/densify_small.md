# Store a Small Design as Base Matrices

Converts the sparse blocks of a design to base matrices when the
equations together carry fewer than `min_dim` coefficients and no block
is sparse because the caller asked for it.

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
records the same crossover for the factorization. Measured on
[`MASS::Cars93`](https://rdrr.io/pkg/MASS/man/Cars93.html), a lasso over
twelve standardized columns beside `random(~ 1 | Manufacturer)` (46
coefficients): the fit at a held \\\lambda\\ takes 2.55 s with the
design dense against 6.69 s sparse, and a path of ten values 30.1 s
against 67.8 s, with the coefficients agreeing to 6e-12 and the same
\\\lambda\\ chosen. The conversion is applied only where no block moves
with the coefficients and no structural term is present, since those
designs are rebuilt during the fit and a rebuilt block would come back
sparse. A block the caller asked to be sparse (`sparse = TRUE` on a
term, `linpar_options(sparse = TRUE)`, a sparse matrix as input) is
never converted, and then neither is the rest of the design: the design
keeps the storage it was given.
