# The Hyperparameters That the Outer Criterion Estimates

Lists the hyperparameters of a fit that
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
can hold, one row for each coordinate of the outer search, with the
address that `which` uses, the estimate on the scale of the
hyperparameter and on the scale of its link.

## Usage

``` r
profile_hyper_table(fit, spec, blocks)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec:

  The specification of the fit.

- blocks:

  Its blocks.

## Value

A data frame with the columns `key`, `parameter`, `term`, `name`,
`estimate`, `eta` and `path`, and the attributes `links` and `members`,
those of the outer index extended by the path's rows. It has no rows
where the fit estimates no hyperparameter.

## Details

The hyperparameters of a kinked penalty that a path selects follow those
of the outer search, with `path` `TRUE`. They have no coordinate in the
outer criterion, so
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
reads them through the grid only.
