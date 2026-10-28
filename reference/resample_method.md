# Resampling Schemes for the Bootstrap of a Fit

Four schemes generate the data sets from which
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
refits a model.
[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md)
draws a new response from the fitted model,
[`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md)
draws observations with replacement,
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
draws groups of observations with replacement, and
[`bayesian_resample()`](https://statmodels7.github.io/statmodels7/reference/bayesian_resample.md)
gives the observations random positive weights. Each is an object, and
the object is passed as the argument `resample` of
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md),
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
and
[`statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.md).

`resample_method` is the abstract parent of the four classes. It cannot
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
variance, which the parametric bootstrap does not. The Bayesian scheme
is a smoothed version of the case scheme: its weights are positive, so
no replicate loses a level of a factor or separates a logistic
regression, and its variance is close to that of the case scheme.

With a kinked penalty
([`modelterms7::lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.html),
[`modelterms7::enet()`](https://statmodels7.github.io/modelterms7/reference/enet.html),
[`modelterms7::scad()`](https://statmodels7.github.io/modelterms7/reference/scad.html),
[`modelterms7::mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.html))
the replicates are those of the penalized estimator. They carry its bias
toward zero, so an interval taken from them covers the expected value of
that estimator and not the coefficient, and the frequency with which a
coordinate is nonzero in the replicates, which
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
prints, is not the probability that the fit selects that coordinate.
Every scheme draws \\n\\ observations; resampling fewer than \\n\\ is
not offered.

A resampled data set is represented by frequency weights, which multiply
the weights of the fit. The design is not rebuilt, so the knots of a
basis and the levels of a factor stay those of the fit, and a refit on
the weights is the refit on the duplicated rows with the basis held. A
cluster scheme beside a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term is the exception: the design is rebuilt on the rows drawn, with the
`build_info` lists of the fit for every term other than the
random-effect terms (see
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)).

## See also

[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md),
[`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md),
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
and
[`bayesian_resample()`](https://statmodels7.github.io/statmodels7/reference/bayesian_resample.md)
for the four schemes,
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
