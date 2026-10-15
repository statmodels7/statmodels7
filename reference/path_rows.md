# Which Hyperparameters a Path Has to Select

The rows of an index like
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)'s,
for the hyperparameters whose penalty has a kink.

## Usage

``` r
path_rows(spec, blocks, hyper)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- blocks:

  The blocks.

- hyper:

  The hyperparameters.

## Value

A data frame with one row per axis, with columns `parameter`, `term` and
`name`, and a list column `members` that holds, for each axis, the data
frame of the hyperparameters it stands for. It has zero rows where no
hyperparameter is to be estimated.

## Details

A gradient search is not the instrument for these. The penalized mode is
a piecewise smooth function of the hyperparameter, differentiable while
the active set holds and turning a corner every time a coefficient joins
it or leaves, so a criterion read at that mode inherits the corners and
a quasi-Newton step is reading a slope that is about to change. A grid
does not care, and warm starts make it cheap.

Every hyperparameter of a kinked penalty is varied, except those that
its term holds and those of an `id` group that one member holds.
