# The Covariance Blocks That a Summary Prints Ahead of the Equations

Returns one block record per covariance class shared by more than one
member: the standard deviations and correlations of the shared prior, in
which each coordinate is named for the equation, the term and the column
to which it belongs.

## Usage

``` r
summary_class_blocks(
  spec,
  design,
  tables,
  edf = NULL,
  coef = NULL,
  hyper = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- tables:

  The per-parameter block lists.

- edf:

  The per-term degrees of freedom, or `NULL`.

- coef:

  The coefficients, or `NULL`. With `hyper`, they let the block report
  what the class itself spends rather than what its members' terms do.

- hyper:

  The hyperparameters, or `NULL`.

## Value

A list of block records in the shape that
[`print_block()`](https://statmodels7.github.io/statmodels7/reference/print_block.md)
reads. Each record carries its coordinates in `coords`. The list is
empty where every class has fewer than two members.

## Details

A shared block belongs to no single member. Printed inside the block of
one member, it would report the standard deviation of an effect on
`sigma` under a term of `mu`, and the block of the other member would
state that there is nothing to report on its own. The block is therefore
printed once, before the equations, where its coordinates can be given
names.

The rows are those that
[`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
held apart, so the numbers are produced in one place; this function adds
the heading of the class and the legend that gives the term of each
coordinate.

## See also

[`class_coords()`](https://statmodels7.github.io/statmodels7/reference/class_coords.md)
for the naming,
[`class_notes()`](https://statmodels7.github.io/statmodels7/reference/class_notes.md)
for the note that gives the same information in prose.
