# Which Rows of a Table a Summary Prints

Every row of a short table, and of a long one the hyperparameters
together with the first few coefficients, the rest reported as a count.

## Usage

``` r
block_rows_shown(tb, max_coef = NULL, cap = 12L, show = 10L)
```

## Arguments

- tb:

  A summary table.

- max_coef:

  How many coefficient rows a long table keeps, `Inf` or `NA` for all of
  them. `NULL` reads the option `statmodels7.summary_max_coef`, and then
  `show`.

- cap:

  The number of rows up to which a table is printed whole. A table is
  cut only above `max(cap, max_coef + number of hyperparameter rows)`.

- show:

  The number of coefficient rows kept when `max_coef` is `NULL` and the
  option is not set.

## Value

An integer vector of row positions.

## Details

A block of a few coefficients is printed whole. A long block is abridged
to its first `max_coef` coefficient rows, and the rows that are dropped
remain in [`coef()`](https://rdrr.io/r/stats/coef.html). The
hyperparameters are never dropped, whatever the length, because they
govern every coefficient under them.
