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

The question is asked of the term's own specification (`spec$linear`)
and never of a suffix in a coefficient's name, a name being a label and
this is a fact about the construction.
