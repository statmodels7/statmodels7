# Resolve the Verbosity Setting

Turns a level or a named logical vector into the four switches `outer`,
`blocks`, `inner` and `optimizer` that the fit reads.

## Usage

``` r
verbosity(verbose)
```

## Arguments

- verbose:

  A number from 0 to 3, or a named logical vector with any of `outer`,
  `blocks`, `inner`, `optimizer`.

## Value

A list of four logicals named `outer`, `blocks`, `inner` and
`optimizer`.

## Details

Each level adds output: `1` shows the outer search and the alternation
(switches `outer` and `blocks`), `2` adds the inner method's own
iterations (`inner`), and `3` adds the traces of the optimizers
(`optimizer`). The named form sets the four switches independently,
which allows the alternation to be shown while a chatty inner optimizer
is silent. A switch that is left out is off.
