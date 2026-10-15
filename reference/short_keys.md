# A Term Key Shortened for Display

The key with the arguments of its leading call reduced to the first one,
so that a trace line names the term instead of repeating its whole
specification.

## Usage

``` r
short_keys(x)
```

## Arguments

- x:

  A character vector of keys.

## Value

A character vector the same length.

## Details

A key is the call that produced the term, and for a penalty over a
sub-term it is that call followed by `::parameter::sub-term`. The call
is what grows when a term has many arguments, so that printing it in
full makes a line of an outer trace unreadable. Only the leading call is
shortened, and only past its first argument, so that
`s(x, bspline_smooth(k = 20))` and `s(z, bspline_smooth(k = 8))` stay
apart. Everything after `::` is kept whole, because that distinguishes
one entry of a term from another.

Where shortening would make two labels the same, the full labels are
returned, all of them, so that every line of a trace follows one
convention.
