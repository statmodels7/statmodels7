# Is a Criterion Minimized?

Reports whether a criterion is to be made small: `TRUE` for
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
which estimate prediction error, and `FALSE` for
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
which are log-likelihoods and are maximized.

## Usage

``` r
outer_minimize(method)
```

## Arguments

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md),
  as
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  or [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  returns it. Only its `kind` is read.

## Value

A single logical.

## Details

The outer search always minimizes, so this decides whether the criterion
goes in as it stands or with its sign turned. Asking the method, rather
than writing the list of kinds into the search, means a criterion added
later declares its own orientation in one place.

## See also

[`outer_k()`](https://statmodels7.github.io/statmodels7/reference/outer_k.md)
for the other property of a criterion the search reads,
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
for the criteria this answers `TRUE` for.
