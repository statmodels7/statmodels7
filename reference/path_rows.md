# Which Hyperparameters a Path Has to Select

The rows of an index like
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)'s,
for the hyperparameters whose penalty has a kink.

## Usage

``` r
path_rows(spec, blocks, hyper, method)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- blocks:

  The blocks.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Value

A data frame of `parameter`, `term` and `name`.

## Details

A gradient search is not the instrument for these. The penalized mode is
a piecewise smooth function of the hyperparameter, differentiable while
the active set holds and turning a corner every time a coefficient joins
it or leaves, so a criterion read at that mode inherits the corners and
a quasi-Newton step is reading a slope that is about to change. A grid
does not care, and warm starts make it cheap.

Which hyperparameters are varied is read from the penalty by
[`kink_hypers()`](https://statmodels7.github.io/statmodels7/reference/kink_hypers.md)
unless the method names them.
