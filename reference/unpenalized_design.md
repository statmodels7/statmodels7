# The Unpenalized Columns of an Equation

Returns the columns of an equation's design that no penalty of their
term reaches.

## Usage

``` r
unpenalized_design(terms, X = NULL, blocks = NULL, n)
```

## Arguments

- terms:

  A named list of built additive terms of one equation.

- X:

  The equation's design, the terms' blocks side by side, or `NULL` to
  take each block from
  [`modelterms7::term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.html).

- blocks:

  The column range of each term in `X`, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  records it; ignored when `X` is `NULL`.

- n:

  The number of rows.

## Value

A dense matrix with `n` rows, possibly with no column.

## Details

Each built additive term reports its penalties through
[`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html),
whose `index` gives the columns of the term's block that a penalty
covers. The remaining columns are free, and they are the ones that can
carry a constant without a penalty identifying it. A random intercept
spans the constant, but its columns are penalized, so it does not count
as an intercept here.

## See also

[`spans_constant()`](https://statmodels7.github.io/statmodels7/reference/spans_constant.md),
[`statmod_held_levels()`](https://statmodels7.github.io/statmodels7/reference/statmod_held_levels.md)
