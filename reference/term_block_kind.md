# What Kind of Block a Term Reports As

Which of the four readings a term gets in a summary: its coefficients, a
smooth's linear part and smoothing parameter, a random effect's variance
parameters, or a selection's survivors.

## Usage

``` r
term_block_kind(term)
```

## Arguments

- term:

  A built term.

## Value

One of `"structural"`, `"breakpoint"`, `"parametric"`, `"smooth"`,
`"random"`, `"selection"`, `"penalized"`.

## Details

The classification is by the term's class and by its penalties, not by
its label, so a term given a name of its own is read the same way. The
penalties are the ones the term declares through
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html),
so a term penalized over part of its parameters – a segmented term's
changes, a filter's deviations – is read as penalized rather than as
parametric, and is a selection when any of its penalties has a kink.
