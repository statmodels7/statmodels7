# The Hyperparameters an Outer Method Estimates

One row per hyperparameter of a penalized term that is fitted in the
joint system, with the link that carries it onto the whole line.

## Usage

``` r
outer_hyper_index(spec, blocks)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- blocks:

  The block split, as
  [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns it.

## Value

A data frame with one row per estimated hyperparameter and three
character columns, `parameter`, `term` and `name`. The links that carry
each onto the whole line ride on the `"link"` attribute, a list the same
length as the frame has rows, since a list column would not survive the
subsetting the search does.

## Details

A term whose penalty has a kink is left out. Its coefficients are
estimated by a coordinate descent with everything else held fixed, and
the criterion is a Laplace approximation, which asks for a second
derivative that does not exist at a kink. Those hyperparameters go to a
path instead.

A hyperparameter its own term holds is left out too, whatever kind of
penalty it belongs to:
[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)
is what says so.

The result is the index the whole outer search runs on. Its row count is
the dimension of that search, and at zero rows the criterion does not
run at all.

## See also

[`hyper_to_eta()`](https://statmodels7.github.io/statmodels7/reference/hyper_to_eta.md)
for the map those links define,
[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)
for what is excluded.
