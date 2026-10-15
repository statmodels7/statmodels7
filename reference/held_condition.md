# The Condition a Held Coordinate Warns Through

A classed warning, so that a caller can silence or catch the one raised
when
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
holds a coordinate about which the information carries nothing, without
silencing every other warning of a fit.

## Usage

``` r
held_condition(msg)
```

## Arguments

- msg:

  The message.

## Value

A condition of class `statmod_held_coord`.

## See also

[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md),
which finds them, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which raises this.
