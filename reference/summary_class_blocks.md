# The Covariance Blocks a Summary Prints Ahead of the Equations

One block record per covariance class spanning more than one term: the
standard deviations and correlations of the shared prior, with each
coordinate named for the equation, the term and the column it belongs
to.

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

A list of block records, in the shape
[`print_block()`](https://statmodels7.github.io/statmodels7/reference/print_block.md)
reads, each carrying its coordinates in `coords`. Empty where no class
spans more than one term.

## Details

A shared block is the property of no single term. Printed inside one
member's block – which is where it was, under whichever member the walk
reached first – it reports the standard deviation of an effect on
`sigma` under a term of `mu`, and the other member's block says there is
nothing to report on its own. Both statements are true of the term and
neither is what a reader wants, so the block is lifted out and printed
once, before the equations, where its coordinates can be given names.

The rows are the ones
[`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
held apart, so the numbers are produced in exactly one place; what is
added here is the class's own heading and the legend saying which term
each coordinate came from.

## See also

[`class_coords()`](https://statmodels7.github.io/statmodels7/reference/class_coords.md)
for the naming,
[`class_notes()`](https://statmodels7.github.io/statmodels7/reference/class_notes.md)
for the note that says the same thing in prose.
