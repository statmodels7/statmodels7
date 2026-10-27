# Evaluate the Outer Criterion at Given Hyperparameters

Returns a function of the free coordinates of the hyperparameters that
refits the model there and returns the criterion and the effective
degrees of freedom of the term.

## Usage

``` r
criterion_evaluator(fit, cd)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- cd:

  The list returned by
  [`criterion_setup()`](https://statmodels7.github.io/statmodels7/reference/criterion_setup.md).

## Value

A function of the free coordinates that returns a list with `value` and
`edf`.

## Details

The refit is the one of
[`criterion_differenced()`](https://statmodels7.github.io/statmodels7/reference/criterion_differenced.md):
the inner fit is warm-started at the coefficients of `fit`, with the
inner method, the budget and the information of the fit, and the
criterion is read at the refitted mode by the function the outer search
uses
([`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md),
or
[`statmod_pe()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe.md)
for
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)).
The value is `NA` where the refit gives no usable criterion.
