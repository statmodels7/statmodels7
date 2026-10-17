# The Value at Which Each Sharing Label Is Held

Returns one entry, keyed by the label, for each sharing label whose
hyperparameter a member holds at a value.

## Usage

``` r
statmod_held_ids(units)
```

## Arguments

- units:

  The penalized units, as
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
  returns them.

## Value

A named list of numbers, empty where no label is held.

## Details

The members of a group estimate one value, so a value written on one of
them is the value of the group. Two members held at different values
contradict each other, and the function signals an error that names both
terms.
