# Which Coefficients of a Smooth Are the Linear Part

`TRUE` for the columns a Demmler-Reinsch smooth carries its linear
effect in, which are the ones worth printing.

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

The rest of the block are coefficients of an orthonormal basis of the
wiggly part; individually they say nothing, and what they say jointly is
the effective degrees of freedom, which the block header reports
instead.

The question is asked of the construction and never of a suffix in a
coefficient's name, a name being a label. What the construction says is
the **penalty**: a column the roughness matrix leaves alone is a column
of the null space the smoother kept, which is what carries the linear
effect. Read that way the answer follows the smoother rather than a
flag: an order-2 penalty leaves one column free, an order-3 penalty two,
and a periodic basis none, its null space being the constant alone and
the constant belonging to the model's intercept.

It replaces a reading of `spec$linear`, which was the term's record of
the same fact while a smooth was always a B-spline with a
second-derivative penalty. Under a factor `by` the block is one copy per
level and only the first level's column is marked, which is what that
reading did too.
