# The Rows of a Block Formatted for Printing

The six numeric columns of a summary table rendered as strings, with the
cells that have no number for a hyperparameter row left empty and the
name of each such row carrying the origin of its value.

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

  What to head the statistic column with.

## Value

A list with `cells`, a character matrix of six columns, and `name`, the
row labels.

## Details

A hyperparameter row prints numbers where there are any: one estimated
by a marginal criterion carries a standard error and an interval. Where
there is none the columns are blank. What put the value there goes in
the name of every hyperparameter row, in square brackets (the role
`"fixed"` or the criterion that chose the value).

The whole block is formatted in one call, its own rows and every
compartment's together, so that the widths to which a column is padded
are the same throughout and the compartments line up with the table
above them.
