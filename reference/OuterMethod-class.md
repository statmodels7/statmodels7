# How the Hyperparameters Are Estimated

The criterion
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
build: which subspace of the coefficients is integrated over, and which
information matrix enters the determinant.

## Usage

``` r
OuterMethod(
  kind = character(0),
  hessian = character(0),
  k = integer(0),
  n_values = integer(0),
  min_ratio = integer(0),
  nfolds = integer(0),
  rule = character(0),
  folds = integer(0),
  over = character(0)
)
```

## Arguments

- kind:

  `"reml"`, `"ml"`, `"aic"`, `"bic"` or `"cv"`.

- hessian:

  `"expected"` or `"observed"`.

- k:

  The price of one degree of freedom, for a prediction-error criterion.
  `NA` where the method resolves it against the sample size.

- n_values:

  How many points a path over a kinked hyperparameter visits.

- min_ratio:

  The smallest kink a path reaches, as a fraction of the one that
  empties the block.

- nfolds:

  How many folds cross-validation uses.

- rule:

  `"min"` or `"1se"`.

- folds:

  A fold number per observation, or `integer(0)`.

- over:

  Which hyperparameters a path varies, or `character(0)` for the ones
  that set the size of the kink.

## Value

An object of class `OuterMethod`.

## See also

[`reml`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`cv`](https://statmodels7.github.io/statmodels7/reference/cv.md),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
reml()
#> <REML>  observed information
ml(hessian = "observed")
#> <ML>  observed information
aic()
#> <AIC>  observed information
cv(nfolds = 5)
#> <CV>  5 folds, 25 values, rule "min"
```
