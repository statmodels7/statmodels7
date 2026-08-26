# Render a Conflict Report

Formats the result of
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
as one line per masked name, in the shape the attach message prints:
`x pkgA::name masks pkgB::name`. The winner is named on the left of
`masks` and every package it hides on the right, comma separated.

## Usage

``` r
format_conflicts(conflicts)
```

## Arguments

- conflicts:

  A named list as
  [`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
  returns: one entry per masked name, each a character vector of at
  least two package names with the winner first. An empty list gives an
  empty result. Not validated; a vector of length one would render as
  masking nothing.

## Value

A character vector with one element per entry of `conflicts`, in the
order they arrived, with no names. `character(0)` for an empty input, so
a caller can [`c()`](https://rdrr.io/r/base/c.html) the result into a
message unconditionally.

## See also

[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
for the input,
[`statmodels7_attach_message()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_attach_message.md)
for the other half of the same message.
