# The Names of the Structural Tail of the Joint Information

Returns the free parameters of every structural term, under the names
that a coefficient of the same term would carry.

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

A level that an intercept in the same equation carries is held and is
not in the joint information, so it is not here either. The tail
consists of the free parameters in the order of the term, which is the
order in which the information was assembled.
