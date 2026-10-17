# Drop the Prefix That a Set of Coefficient Names Share

Removes the first dotted piece where every name carries the same one and
removing it leaves every name non-empty.

## Usage

``` r
drop_common_prefix(nms)
```

## Arguments

- nms:

  The names.

## Value

The names, shortened where they share a prefix.

## Details

A term composes its coefficients' names from its own and its
parameters', so inside the block of one term the leading piece repeats
on every row and repeats what the heading has already stated. It is
dropped only for printing; [`coef()`](https://rdrr.io/r/stats/coef.html)
and the summary's own tables keep the names with which the fit was
built, which are the names that another call can use as an index.

Only one piece is dropped, and never every piece that the names share:
the term's own name is one, and a set of coefficients that agree further
along, such as `r.1` and `r.2` of a matrix column, would otherwise be
left as bare numbers.
