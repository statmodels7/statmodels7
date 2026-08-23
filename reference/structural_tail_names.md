# The Names of the Structural Tail of the Joint Information

The free parameters of every structural term, under the names a
coefficient of the same term would carry.

## Usage

``` r
structural_tail_names(spec, design)
```

## Arguments

- spec:

  The fitted specification.

- design:

  The design.

## Value

A character vector, empty where there is no structural term.

## Details

A level an intercept in the same equation carries is held and is not in
the joint information, so it is not here either; the tail is the free
ones in the term's own order, which is the order the information was
assembled in.
