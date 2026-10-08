# The Predicted Random Effects of a Term, Summarized

The six-number summary of
[`base::summary()`](https://rdrr.io/r/base/summary.html) (minimum,
quartiles, mean and maximum) of a random-effect term's predicted
effects, one row per coordinate of the effect: the intercept, and each
covariate of a random slope.

## Usage

``` r
effect_spread(term, est)
```

## Arguments

- term:

  A built random-effect term.

- est:

  Its predicted effects, in the order of its coefficients.

## Value

A list with `table`, a data frame with one row per coordinate and the
six summary columns, `levels`, the number of groups, and `group`, the
grouping expression; `NULL` where the term carries no grouping.

## Details

The coefficients of one group are adjacent in the block, so the
predictions are read as a matrix with one row per level and one column
per coordinate, the coordinates named as
[`modelterms7::term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.html)
names them.
