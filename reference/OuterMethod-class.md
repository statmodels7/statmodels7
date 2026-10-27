# How the Hyperparameters Are Estimated

Every outer criterion is an instance of this class:
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
all return one. It carries which subspace of the coefficients is
integrated over, which information matrix enters the determinant, and
the settings that a path or a cross-validation needs.

One class serves every criterion, so a property that a given criterion
does not read is present and unused. The class validates the values of
its properties. The constructors
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
also check `information`, match `rule` against its choices and check `k`
and `marginal`.

## Usage

``` r
OuterMethod(
  kind = character(0),
  hessian = distributions7::Information(),
  k = integer(0),
  nfolds = integer(0),
  rule = character(0),
  folds = integer(0),
  marginal = NULL
)
```

## Arguments

- kind:

  `"reml"`, `"ml"`, `"aic"`, `"bic"` or `"cv"`. Determines which
  criterion is evaluated and, through
  [`outer_minimize()`](https://statmodels7.github.io/statmodels7/reference/outer_minimize.md),
  whether it is made small or large.

- hessian:

  An object that inherits from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html):
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html),
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  or
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md).
  Which information enters the Laplace determinant.
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  settles
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  against the family through
  [`outer_resolve()`](https://statmodels7.github.io/statmodels7/reference/outer_resolve.md).
  The constructors set it from their argument `information`.

- k:

  The price of one effective degree of freedom, for a prediction-error
  criterion. It is `NA_real_` for
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  which resolves it to \\\log n\\ at fit time, and for the criteria that
  do not charge per degree of freedom.

- nfolds:

  How many folds
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  uses. Only `"cv"` reads it.

- rule:

  `"min"` or `"1se"`. Only `"cv"` reads it.

- folds:

  A fold number per observation, or `integer(0)` for folds drawn at fit
  time. Only `"cv"` reads it.

- marginal:

  Which distribution parameters have their unpenalized coefficients
  estimated on the marginal criterion: `NULL` for every parameter except
  the position, `"none"`, `"all"`, or a vector of names. Only `"reml"`
  and `"ml"` read it; see
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).

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
ml(information = observed_info())
#> <ML>  observed information
aic()
#> <AIC>  observed information
cv(nfolds = 5)
#> <CV>  5 folds, rule "min"
```
