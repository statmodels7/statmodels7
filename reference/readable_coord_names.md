# A Coordinate's Name in Place of Its Number

Rewrites the names a multivariate family gives the quantities of a
matrix parameter – `sd_v1`, `cor_v1_v2` – into the same quantities said
of the coordinates this model is written in: `sd[mu:(Intercept)]`,
`cor[mu:(Intercept), sigma:(Intercept)]`.

## Usage

``` r
readable_coord_names(nm, labels)
```

## Arguments

- nm:

  The names as the family gives them.

- labels:

  One label per coordinate, in the matrix's own order.

## Value

A character vector as long as `nm`.

## Details

The family names its coordinates by position because that is all it has,
and every one of its readings – a standard deviation, a correlation, a
partial correlation, a conditional scale – is named `prefix_vi` or
`prefix_vi_vj`. The rewrite is on that shape alone: a name it does not
match, or one whose index no label answers for, is left exactly as it
is, so a family declaring a reading of another kind keeps its own name
rather than being renamed into a wrong one.

## See also

[`class_coords()`](https://statmodels7.github.io/statmodels7/reference/class_coords.md)
and
[`term_coord_labels()`](https://statmodels7.github.io/statmodels7/reference/term_coord_labels.md)
for the labels,
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for the names being rewritten.
