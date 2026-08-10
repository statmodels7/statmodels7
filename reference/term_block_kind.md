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

One of `"parametric"`, `"smooth"`, `"random"`, `"selection"`,
`"penalized"`.

## Details

The classification is by the term's class and by its penalty, not by its
label, so a term given a name of its own is read the same way.
