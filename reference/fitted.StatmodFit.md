# The Fitted Values of a Model

One distribution parameter's fitted values, as a vector with one value
per row of the fitting data.

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

  Not used. Any further argument signals an error, since the values at
  new data are given by `predict(fit, what, newdata = new_data)`.

## Value

A numeric vector with one value per row of the fitting data, that
parameter's fitted values on its own scale.

## Details

The result is one vector, never the whole set, so that this and
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md)
are a matched pair and a diagnostic drawn from them needs no unpacking.

The default is the first parameter, not the mean. A family may have no
mean, a Cauchy being one, whereas every family has a first parameter.

The whole set at once is `predict(fit, "parameter")`, and the mean,
where it exists, is `predict(fit, "mean")`.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for the other quantities and for new data,
[`residuals.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/residuals.StatmodFit.md)
for the matched diagnostic.
