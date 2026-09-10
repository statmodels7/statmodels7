# What the Search Reports About Itself

The one place the search's own verdict is worded, read by both
[`print.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodFit.md)
and
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
so the two cannot drift apart.

It says whether the outer search met its stopping rule, which is a
property of the SEARCH and not of the point it stopped at. The point is
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)'s
question, and it is the one a reader has.

## Usage

``` r
search_verdict(converged)
```

## Arguments

- converged:

  The optimizer's own flag.

## Value

A single string.
