# Where a Penalized Unit's Coordinates Are

The positions the unit's penalty covers, in whichever vector it is
addressed in: the stacked coefficients for an ordinary unit, the
structural term's own parameters for one whose coefficients are a
filter's, and the joint vector of coefficients and free structural
parameters for a mixed class.

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
or which of a term's components they belong to, needs the same answer
for every kind of unit, and the kinds keep their positions in different
fields. An ordinary unit and a covariance class among the coefficients
carry them in `index`. A structural unit has no `index` and carries them
in `cols`, among the term's own parameters. A mixed class has neither
and carries them in `joint`.

## See also

[`unit_beta()`](https://statmodels7.github.io/statmodels7/reference/unit_beta.md)
for the values at those positions.
