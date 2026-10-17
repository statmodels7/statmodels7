# A Coordinate's Name in Place of Its Number

Rewrites the names that a multivariate family gives to the quantities of
a matrix parameter (`sd_v1`, `cor_v1_v2`) into names for the same
quantities in terms of the coordinates in which this model is written:
`sd[mu:(Intercept)]`, `cor[mu:(Intercept), sigma:(Intercept)]`.

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

A multivariate family names its coordinates by position, and each of its
readings – a standard deviation, a correlation, a partial correlation, a
conditional scale – is named `prefix_vi` or `prefix_vi_vj`. The rewrite
applies only to that shape: a name that does not match, or whose index
has no label, is left unchanged, so a family declaring a reading of
another kind keeps its own name.

## See also

[`class_coords()`](https://statmodels7.github.io/statmodels7/reference/class_coords.md)
and
[`term_coord_labels()`](https://statmodels7.github.io/statmodels7/reference/term_coord_labels.md)
for the labels,
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for the names that are rewritten.
