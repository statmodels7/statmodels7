# The Sub-Term to Which a Penalty Entry Belongs

Returns the term whose own columns an entry of
[`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html)
covers: the term itself where the entry is over its own block, and the
sub-term developing one of its parameters where it is not.

## Usage

``` r
entry_owner(term, ii)
```

## Arguments

- term:

  One built term.

- ii:

  The entry's columns, in the term's own block.

## Value

One built term.

## Details

A term that develops a parameter over another term reports that term's
penalties as its own, and each entry names the columns that it covers.
Those columns are the sub-term's block, so a question about the penalty
(above all, what its coordinates are) is a question about the sub-term,
and reading it from the parent returns nothing: `nl()` has no grouping,
and a correlated random effect inside it does.

The walk recurses, translating the entry's positions into the sub-term's
own on the way down, so a development two levels deep is reached.

## See also

[`term_coord_labels()`](https://statmodels7.github.io/statmodels7/reference/term_coord_labels.md),
its caller.
