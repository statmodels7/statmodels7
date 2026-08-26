# The Fitted Values of a Model

One distribution parameter's fitted values, as a vector of the data's
length.

## Usage

``` r
# S3 method for class 'StatmodFit'
fitted(object, what = NULL, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- what:

  Which distribution parameter, a string naming one of the family's, or
  `NULL` for the first.

- ...:

  Unused.

## Value

A numeric vector of length `nobs(object)`, that parameter's fitted
values on its own scale.

## Details

The result is one vector, never the whole set, so that this and
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md)
are a matched pair and a diagnostic drawn from them needs no unpacking.

The default is the **first** parameter, never the mean. A family may
have no mean, a Cauchy being one, and a default that fails on a
legitimate family is worse than one that always answers.

The whole set at once is `predict(fit, "parameter")`, and the mean,
where it exists, is `predict(fit, "mean")`.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for the other quantities and for new data,
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md)
for the matched diagnostic.
