# Data-Based Intercepts of Families That Supply Them

Overwrites, in the list that
[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
builds, the intercept of each parameter for which
[`distributions7::distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.html)
returns a value, carried onto that parameter's own link.

## Usage

``` r
intercept_start_from_data(spec, eta)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- eta:

  The named list of link-scale intercepts.

## Value

`eta`, in which the answered entries are replaced.

## Details

The intercept-only fit is the right start for almost every family. It is
the wrong one for the mixing weight of
[`distributions7::zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.html),
and for `skewnormal1_distrib()` and the pseudo-Huber families, for which
[`distributions7::distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.html)
has methods of their own. For the mixing weight, the overdispersion of
the parent can absorb the excess zeros without covariates, the weight
goes to the edge of its domain, and the link is flat there, so that a
model with covariates started from it does not leave the edge. The
family returns the observed proportion of zeros, and the value is mapped
through the parameter's link, so a probit or a cloglog chart has its own
linear predictor and not a logit.

## See also

[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
