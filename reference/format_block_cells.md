# The Rows of a Block Formatted for Printing

The six numeric columns of a summary table rendered as strings, with the
cells a hyperparameter row has no number for left empty and the name of
each such row carrying what put the value there.

## Usage

``` r
format_block_cells(tb, digits = 4L)
```

## Arguments

- tb:

  A summary table.

- digits:

  Significant digits.

## Value

A list with `cells`, a character matrix of six columns, and `name`, the
row labels.

## Details

A hyperparameter row prints numbers where there are any: one estimated
by a marginal criterion carries a standard error and an interval. Where
there is none the columns are blank. What put the value there goes in
the NAME, on every hyperparameter row rather than only on the ones with
nothing else in them: written into the column where a standard error
would have been it marked a held or path-chosen row and never a REML
one, whose column is occupied, so the note at the foot spoke of a mark
that was never printed.

The whole block is formatted in ONE call, its own rows and every
compartment's together, so that the widths a column is padded to are the
same throughout and the compartments line up under the table they sit
beneath.
