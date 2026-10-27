# The Failures Recorded in a Bootstrap Entry

Counts, by cause, the draws that failed before the last replicate kept.

## Usage

``` r
boot_failures(x, n = NULL)
```

## Arguments

- x:

  A `StatmodBootstrap` entry.

- n:

  Count only the draws up to the `n`-th replicate kept; `NULL` counts
  all of them.

## Value

An integer vector named `error`, `mode_error` and `unidentified`.
