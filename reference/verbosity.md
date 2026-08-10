# Resolve the Verbosity Setting

Turns a level or a named logical vector into the three switches the fit
reads.

## Usage

``` r
verbosity(verbose)
```

## Arguments

- verbose:

  A number from 0 to 3, or a named logical vector with any of `blocks`,
  `inner`, `optimizer`.

## Value

A list of three logicals.

## Details

The levels name the loops rather than counting them: `1` shows the
alternation, `2` the inner method's own iterations, `3` the optimizers'
traces. The named form exists because the three are genuinely
independent – watching the alternation while silencing a chatty inner
optimizer is the common case.
