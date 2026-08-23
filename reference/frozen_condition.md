# The Condition a Frozen Block Raises

A warning of its own class, so that a caller reporting the same thing in
its own channel can muffle it rather than repeat it.

## Usage

``` r
frozen_condition(msg)
```

## Arguments

- msg:

  The message.

## Value

A condition.

## Details

[`summary.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
calls [`vcov`](https://rdrr.io/r/stats/vcov.html) more than once and
would raise the warning once per call; it says it once, as a note.
