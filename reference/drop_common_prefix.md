# Drop the Prefix a Set of Coefficient Names Share

The first dotted piece, where every name carries the same one and
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
on every row and says what the heading has already said. It is dropped
for the printing alone; [`coef()`](https://rdrr.io/r/stats/coef.html)
and the summary's own tables keep the names the fit was built with,
which are the ones another call can be indexed by.

**one piece** and never every piece they share: the term's own name is
one, and a set of coefficients that happen to agree further along, `r.1`
and `r.2` of a matrix column, would otherwise be left as the bare
numbers.
