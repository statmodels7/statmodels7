# Print One Block of a Model Summary

The heading of a block, the term read at a glance where it is written in
parameters of its own, its own coefficients, and one indented
compartment per parameter developed over covariates.

## Usage

``` r
print_block(b, digits = 4L, max_coef = NULL, stat = "z")
```

## Arguments

- b:

  A block record from
  [`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md).

- digits:

  Significant digits.

- max_coef:

  How many coefficient rows a long block keeps, as
  [`block_rows_shown()`](https://statmodels7.github.io/statmodels7/reference/block_rows_shown.md)
  takes it.

- stat:

  What to head the statistic column with: `"z"` for the estimate over
  its standard error, `"r"` for the signed root of a restricted test.

## Value

`NULL`, invisibly. Called for the printing.

## Details

A term that develops one of its own parameters carries columns that mean
different things: a break-point's population value and its per-group
deviations are not comparable quantities, and a table that stacks them
reads as a list of numbers and no longer as a model. Each developed
parameter is therefore printed as a compartment of its own, headed by
what develops it, with its fixed effects first and then each penalized
sub-term, rendered the way a block of that kind is rendered at the top
level. A random sub-term is headed by its call and its number of
coefficients and reports the hyperparameters of its prior under the
names of an ordinary random() block, followed by the summary of its
predicted effects
([`effect_spread()`](https://statmodels7.github.io/statmodels7/reference/effect_spread.md));
the predictions themselves are in
[`coef()`](https://rdrr.io/r/stats/coef.html).
