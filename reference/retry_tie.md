# Whether Two Fits Reach the Same Point

Returns `TRUE` where the criteria of two fits, or their objectives where
either has no criterion, differ by at most \\10^{-8}(1 + \lvert b
\rvert)\\, where \\b\\ is the second.

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
uses it after a second fit of a score-driven model from the term's own
start: where the two reach the same point, the converged one is kept,
and a difference at the level of rounding does not decide the choice.
