# Which Rows of a Table a Summary Prints

Every row of a short table, and of a long one the hyperparameters
together with the first few coefficients, the rest reported as a count.

## Usage

``` r
block_rows_shown(tb, cap = 12L, show = 10L)
```

## Arguments

- tb:

  A summary table.

- cap:

  The length above which a table is cut.

- show:

  How many coefficient rows a cut table keeps.

## Value

An integer vector of row positions.

## Details

A block of a few coefficients is printed whole: a threshold that cut it
would hide the very numbers a reader opened the summary for. A block of
many is a column of numbers nobody reads to the end, and what is dropped
is still in [`coef`](https://rdrr.io/r/stats/coef.html). The
hyperparameters are never dropped, whatever the length: they govern
every coefficient under them.
