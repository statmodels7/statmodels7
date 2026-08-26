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
is what grows: printed in full, one line of an outer trace carried the
deparsed
`gas(p = 1, q = 1, time = t, by = ~ridge(~id), links = list(...))` three
times over, which is a line no reader can use. Only the leading call is
shortened, and only past its first argument, so `s(x, k = 20)` and
`s(z, k = 8)` stay apart; everything after `::` is kept whole, that
being what distinguishes one entry of a term from another.

Where shortening would make two labels the same the full ones are
returned, all of them: a shorter label that is ambiguous is worse than a
long one, and deciding per label would leave a reader unable to tell
which convention a given line follows.
