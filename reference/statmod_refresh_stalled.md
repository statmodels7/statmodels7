# Have the Refreshable Terms Run Out of Step Control?

`TRUE` when at least one refreshable term reports, through
[`modelterms7::term_stalled()`](https://statmodels7.github.io/modelterms7/reference/term_stalled.html),
that it has not settled and has no step control left, and every other
term asked has either settled or reports the same.
[`fit_working()`](https://statmodels7.github.io/statmodels7/reference/fit_working.md)
reads it to end a working phase that cannot settle.

## Usage

``` r
statmod_refresh_stalled(spec, design, which = "all")
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, whose refresh state holds the terms asked.

- which:

  Which entries to ask, as in
  [`statmod_refresh_settled()`](https://statmodels7.github.io/statmodels7/reference/statmod_refresh_settled.md).

## Value

A single logical, `FALSE` when there is nothing to ask.

## See also

[`statmod_refresh_settled()`](https://statmodels7.github.io/statmodels7/reference/statmod_refresh_settled.md),
[`fit_working()`](https://statmodels7.github.io/statmodels7/reference/fit_working.md)
