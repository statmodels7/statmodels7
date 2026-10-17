# Which Coefficients of a Smooth Are the Linear Part

Returns `TRUE` for the columns of a smooth's block that are not touched
by any penalty of the term, which carry its unpenalized (linear or
polynomial) effect and are the columns that the summary prints.

## Usage

``` r
smooth_linear_cols(term, k)
```

## Arguments

- term:

  A built smooth term.

- k:

  The number of columns in its block.

## Value

A logical vector of length `k`.

## Details

The remaining coefficients of the block are coordinates of an
orthonormal basis of the wiggly part. A single one of them has no
interpretation, and together they determine the effective degrees of
freedom, which the block header reports instead.

The columns are identified from the construction and not from a suffix
in a coefficient's name, which is only a label. The decisive quantity is
the penalty: a column that the roughness matrix leaves alone belongs to
the null space that the smoother kept, and that null space carries the
unpenalized effect. The answer therefore follows the smoother: an
order-2 penalty leaves one column free, an order-3 penalty two, and a
periodic basis none, because its null space contains only the constant
and the constant belongs to the model's intercept. Under a factor `by`
the block is one copy per level and the free columns of every level are
marked, whether the levels share a smoothing parameter or not.
