# Replace an Intercept Where the Family Reads Its Start Off the Data

Overwrites, in the list
[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
builds, the intercept of each parameter
[`distributions7::distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.html)
answers for, carried onto that parameter's own link.

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

`eta`, with the answered entries replaced.

## Details

The intercept-only fit is the right start for almost every family, and
the wrong one for the mixing weight of
[`distributions7::zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.html):
without covariates the parent's overdispersion can absorb the excess
zeros, the weight goes to the edge of its domain, and the link is flat
there, so a model with covariates started from it never leaves –
measured on a negative binomial with a true weight of 0.25, the fit
reported convergence at a weight of zero, 7.6 log-likelihood units below
the interior maximum. The family answers with the observed proportion of
zeros, and the value is mapped through the parameter's link, so a probit
or a cloglog chart gets its own linear predictor and not a logit.

## See also

[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
