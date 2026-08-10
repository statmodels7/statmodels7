# How the Hyperparameters Are Estimated

The criterion
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
build: which subspace of the coefficients is integrated over, and which
information matrix enters the determinant.

## Usage

``` r
OuterMethod(kind = character(0), hessian = character(0))
```

## Arguments

- kind:

  `"reml"` or `"ml"`.

- hessian:

  `"expected"` or `"observed"`.

## Value

An object of class `OuterMethod`.

## See also

[`reml`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
reml()
#> <REML>  expected information
ml(hessian = "observed")
#> <ML>  observed information
```
