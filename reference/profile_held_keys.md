# The Hyperparameters That the Terms Hold

Lists the hyperparameters that the call of a term fixes, and those whose
shared label another member of the group fixes. The fit does not
estimate either, so neither has a profile.

## Usage

``` r
profile_held_keys(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A data frame with the columns `key`, in the format of `which`, and
`why`, a clause that completes the sentence "the key is a hyperparameter
that ...".
