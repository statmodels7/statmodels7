# The Table of Bootstrap Entries

`boot_table()` returns the data frame that
[`statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.md)
prints, and `boot_table_text()` the same as text for an error message.

## Usage

``` r
boot_table(ents)

boot_table_text(ents)

boot_ambiguous(ents)
```

## Arguments

- ents:

  A list of entries.

## Value

A data frame, or a string.
