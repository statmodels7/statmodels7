# Where a Penalized Unit's Coordinates Are

The positions the unit's penalty covers, in whichever vector it is
addressed in: the stacked coefficients for an ordinary unit, and the
structural term's own parameters for one whose coefficients are a
filter's.

## Usage

``` r
unit_positions(u)
```

## Arguments

- u:

  One unit, from
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md),
  or a class carrying the same two fields.

## Value

An integer vector.

## Details

A reader that only needs to know how many coordinates a penalty covers,
or which of a term's components they belong to, wants the same answer
for both kinds and cannot get it from one field: an ordinary unit leaves
`cols` empty where it spans several equations, and a structural one has
no `index` at all.

## See also

[`unit_beta()`](https://statmodels7.github.io/statmodels7/reference/unit_beta.md)
for the values at those positions.
