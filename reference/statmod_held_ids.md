# The Value Each Sharing Label Is Held At

One entry per label a member holds its hyperparameter at, keyed by the
label.

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
them is the group's. Two members held at DIFFERENT values are a
contradiction rather than a preference between them, and are rejected
here, at the one place both are visible.
