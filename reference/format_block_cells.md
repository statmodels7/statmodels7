# The Rows of a Block Formatted for Printing

Renders the six numeric columns of a summary table as strings. The cells
that have no number for a hyperparameter row are left empty, and the
name of each such row carries the origin of its value.

## Usage

``` r
format_block_cells(tb, digits = 4L, stat = "z")
```

## Arguments

- tb:

  A summary table.

- digits:

  Significant digits.

- stat:

  The label at the head of the statistic column.

## Value

A list with `cells`, a character matrix of six columns, and `name`, the
row labels.

## Details

A hyperparameter row prints numbers where there are any: one estimated
by a marginal criterion carries a standard error and an interval. Where
there is none the columns are blank. The name of every hyperparameter
row carries, in square brackets, the origin of its value (the role
`"fixed"` or the criterion that chose the value).

The whole block is formatted in one call, its own rows and every
compartment's together, so that the widths to which a column is padded
are the same throughout and the compartments line up with the table
above them.
