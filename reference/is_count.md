# A Positive Whole Number

Tests that a value is a single finite number, at least one, with no
fractional part, before it is converted with
[`as.integer()`](https://rdrr.io/r/base/integer.html), which would
truncate it.

## Usage

``` r
is_count(v)
```

## Arguments

- v:

  The value.

## Value

`TRUE` or `FALSE`.
