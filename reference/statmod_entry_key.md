# The Key of One of a Term's Penalties

Returns the key under which a hyperparameter row is filed: the term's
name where the term carries one penalty or the entry has no name, and
`term::entry` where the term carries several penalties and the entry is
named.

## Usage

``` r
statmod_entry_key(term, entries, entry)
```

## Arguments

- term:

  The term's name in the formula.

- entries:

  The term's entries, as
  [`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html)
  returns them.

- entry:

  One of them.

## Value

A single string.

## Details

The composition is written once so that every caller reading a
hyperparameter by key uses the same key. A term carrying one penalty is
keyed by its name without the entry's name, whatever the entry is
called.
