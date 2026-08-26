# The Key of One of a Term's Penalties

What a hyperparameter row is filed under: the term's name where the term
carries one penalty, and `term::entry` where it carries several.

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

The composition is written once because two callers reading a
hyperparameter by a key they each compose would agree only by accident.
A term carrying one penalty keys as it always did, so a formula of
ordinary terms is unaffected by the entry names the terms now supply.
