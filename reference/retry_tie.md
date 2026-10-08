# Whether Two Fits Reach the Same Point

`TRUE` where the criteria of two fits, or their objectives where either
has no criterion, agree to a relative \\10^{-8}\\.

## Usage

``` r
retry_tie(a, b)
```

## Arguments

- a, b:

  Two results of the inner fit, each a list with `res` and `crit`.

## Value

`TRUE` or `FALSE`.

## Details

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
reads it after a second fit of a score-driven model from the term's own
start: where the two reach the same point, the converged one is kept, a
difference at the level of rounding deciding nothing.
