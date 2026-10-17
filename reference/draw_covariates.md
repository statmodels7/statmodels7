# Draw One Replicate's Covariates

Evaluates each generator at the observation count and writes the columns
into the data frame.

## Usage

``` r
draw_covariates(covariates, n, data)
```

## Arguments

- covariates:

  The generators.

- n:

  The observation count.

- data:

  The frame into which the columns are written.

## Value

A data frame.

## Details

A generator that returns the wrong number of values is reported and not
recycled, for the reason given for a coefficient function: R would
recycle the vector without a message and the replicate would be of
another model. A column of `data` under the same name is overwritten,
because naming the column in `covariates` requests that it be drawn.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
