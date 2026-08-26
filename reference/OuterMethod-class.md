# How the Hyperparameters Are Estimated

The class every outer criterion is an instance of:
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
all return one. It carries which subspace of the coefficients is
integrated over, which information matrix enters the determinant, and
the settings a path or a cross-validation needs.

One class serves all five, so a property a given criterion does not read
is present and unused. Use the five constructors, which validate; this
raw one does not.

## Usage

``` r
OuterMethod(
  kind = character(0),
  hessian = character(0),
  k = integer(0),
  nfolds = integer(0),
  rule = character(0),
  folds = integer(0)
)
```

## Arguments

- kind:

  `"reml"`, `"ml"`, `"aic"`, `"bic"` or `"cv"`. Decides which criterion
  is evaluated and, through
  [`outer_minimize()`](https://statmodels7.github.io/statmodels7/reference/outer_minimize.md),
  whether it is made small or large.

- hessian:

  `"expected"` or `"observed"`. Which information enters the Laplace
  determinant, and whether the exact gradient is available.

- k:

  The price of one effective degree of freedom, for a prediction-error
  criterion. `NA_real_` for
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  which resolves it to \\\log n\\ at fit time, and for the criteria that
  do not charge per degree of freedom.

- nfolds:

  How many folds
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  uses. Read by `"cv"` alone.

- rule:

  `"min"` or `"1se"`. Read by `"cv"` alone.

- folds:

  A fold number per observation, or `integer(0)` for folds drawn at fit
  time. Read by `"cv"` alone.

## Value

An object of class `OuterMethod` with one property per argument above.

## See also

[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
for the marginal criteria,
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
for the prediction-error ones,
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) for
cross-validation,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for where one is passed.

## Examples

``` r
reml()
#> <REML>  observed information
ml(hessian = "observed")
#> <ML>  observed information
aic()
#> <AIC>  observed information
cv(nfolds = 5)
#> <CV>  5 folds, rule "min"
```
