# The Multisets of Predictor Indices and Their Sub-Multisets

Returns, for \\P\\ predictors and the orders 0 to `order`, the multisets
of each size, the table of sub-multisets that the compiled functions of
[`joint_scale_derivs()`](https://statmodels7.github.io/statmodels7/reference/joint_scale_derivs.md)
read, and the position of each ordered tuple's multiset.

## Usage

``` r
joint_scale_tables(P, order)
```

## Arguments

- P:

  The number of predictors.

- order:

  The highest size.

## Value

A list with `ms` (the multisets, one matrix per size from 0), `tab` (one
integer matrix per order from 1) and `map` (one integer vector per order
from 1).

## Details

The multisets of size \\k\\ are the sorted tuples over \\1, \dots, P\\,
in lexicographic order. The table of order \\m\\ has one row per
multiset \\B\\ and one column per subset of its positions, a bit mask
with position \\j\\ as bit \\j - 1\\; the entry is the zero-based index
of the sub-multiset among the multisets of its size. The map of order
\\m\\ gives, for each tuple in the order of
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) (the first
index running fastest), the one-based index of its multiset. The result
is cached by \\P\\ and `order`.
