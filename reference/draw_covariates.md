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

  The frame to write into.

## Value

A data frame.

## Details

A generator that answers with the wrong length is reported rather than
recycled, for the reason a coefficient function is: R would recycle it
without a word and the replicate would be of another model. A column of
`data` under the same name is overwritten, as a caller asking for that
column to be drawn means.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
