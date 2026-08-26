# The Price of One Degree of Freedom

Returns \\\kappa\\, the price a prediction-error criterion charges for
one effective degree of freedom, resolved against the sample size where
the method left it open.
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
stores `NA` and gets \\\log n\\ here; every other criterion returns the
`k` it was constructed with.

## Usage

``` r
outer_k(method, n)
```

## Arguments

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).
  Its `kind` and `k` are read.

- n:

  The number of observations, a single positive number.

## Value

A single number: \\\log n\\ for
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
and `method@k` for anything else, which is `NA_real_` for
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
none of which charges per degree of freedom.

## Details

The resolution happens at fit time because
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
takes no data and cannot know \\n\\ at construction. A criterion object
is a specification, reusable across models of different sizes.

## See also

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
for where `k` is set,
[`outer_minimize()`](https://statmodels7.github.io/statmodels7/reference/outer_minimize.md)
for the other property read at the same point.
