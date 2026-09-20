# The Rank of the Active Columns

\\\mathrm{rank}(X_A)\\, the degrees of freedom of a fit whose penalty is
flat on the active block (Tibshirani and Taylor, 2012), summed over the
equations the active coordinates belong to.

## Usage

``` r
active_rank(design, act)
```

## Arguments

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- act:

  The active coordinates, as positions in the stacked vector.

## Value

A single number.

## Details

The information is block diagonal over the equations wherever the
family's per-observation information is, and in any case a coordinate is
identified only if its own equation's active columns identify it, so the
rank is read per equation and added. Nothing of size \\n \times \|A\|\\
is formed twice: the columns are subset once and decomposed once.

## See also

[`design_count_exact()`](https://statmodels7.github.io/statmodels7/reference/design_count_exact.md),
whose certificate makes this unnecessary, and
[`block_column_rank()`](https://statmodels7.github.io/statmodels7/reference/block_column_rank.md).
