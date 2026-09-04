# Print One Block of a Model Summary

The heading of a block, the term read at a glance where it is written in
parameters of its own, its own coefficients, and one indented
compartment per parameter developed over covariates.

## Usage

``` r
print_block(b, digits = 4L, n = NULL)
```

## Arguments

- b:

  A block record from
  [`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md).

- digits:

  Significant digits.

## Value

`NULL`, invisibly. Called for the printing.

## Details

A term that develops one of its own parameters carries columns that mean
different things: a break-point's population value and its per-group
deviations are not comparable quantities, and a table that stacks them
reads as a list of numbers and no longer as a model. Each developed
parameter is therefore printed as a compartment of its own, headed by
what develops it, opening with its hyperparameter under a name that says
what the hyperparameter is, and rendering each sub-term the way a block
of that kind is rendered at the top level. A random development reports
the scale of its effects and one line saying how many predictions there
are and how far they spread, the predictions themselves being in
[`coef()`](https://rdrr.io/r/stats/coef.html).
