# Resampling Schemes for the Bootstrap of a Fit

Three schemes generate the data sets from which
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
refits a model.
[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md)
draws a new response from the fitted model,
[`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md)
draws observations with replacement, and
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
draws groups of observations with replacement. Each is an object, and
the object is passed as the argument `resample` of
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md),
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
and
[`statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.md).

`resample_method` is the abstract parent of the three classes. It cannot
be instantiated, and serves to test membership with
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html).

## Usage

``` r
resample_method(name = character(0))
```

## Arguments

- name:

  A single string, the short name of the scheme.

## Value

Nothing: the class is abstract and its constructor signals an error.
Every object that inherits from it has the property `name`.

## Details

The schemes answer two different questions. The parametric scheme
simulates from the fitted model, so its replicates have the variance
that the model attributes to the estimator, and the bootstrap variance
converges to the model-based one (to
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
where the penalties and hyperparameters are held). The case and cluster
schemes resample the data, so they do not depend on the variance
specification of the model: with a heteroscedastic response fitted by a
homoscedastic model the case bootstrap converges to the sandwich
variance, which the parametric bootstrap does not.

A resampled data set is represented by frequency weights, which multiply
the weights of the fit. The design is not rebuilt, so the knots of a
basis and the levels of a factor stay those of the fit, and a refit on
the weights is the refit on the duplicated rows with the basis held.

## See also

[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md),
[`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md)
and
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
for the three schemes,
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
for the function that takes one.

## Examples

``` r
S7::S7_inherits(case_resample(), resample_method)
#> [1] TRUE
try(resample_method(name = "mine"))
#> Error in S7::new_object(S7::S7_object(), name = name) : 
#>   Can't construct an object from abstract class <resample_method>
```
