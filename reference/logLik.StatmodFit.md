# The Log-Likelihood of a Fit

Returns the value that R's convention expects, with the degrees of
freedom and the number of observations as attributes.
[`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
is a different function: it evaluates the model at parameters and data
chosen by the caller.

## Usage

``` r
# S3 method for class 'StatmodFit'
logLik(object, loglik_type = c("conditional", "marginal"), ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- loglik_type:

  `"conditional"` (the default) or `"marginal"`.

- ...:

  Unused.

## Value

A `logLik` object: a single number with attributes `df`, the degrees of
freedom, and `nobs`, the sum of the prior weights (the number of rows
where the weights are all 1).
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) read it, so
`AIC(fit)` is the conditional AIC.

## Which likelihood, and why it matters

The default is the conditional one: the log-density at the fitted
coefficients, including the penalized coefficients, paired with the
effective degrees of freedom \\\mathrm{tr}\[(H+S)^{-1}H\]\\. A criterion
built on that pair is the conditional AIC, which measures how well the
model describes the groups, curves and states actually observed.

A mixed-model package reports the marginal likelihood instead: the
random effects are integrated out and the count is the number of
estimated parameters, including the variance components. That likelihood
concerns the population, and the AIC built on it is a different number,
not comparable with the conditional one.

Each likelihood keeps its own count. A marginal likelihood against an
effective count, or a conditional one against a parameter count, mixes
the two conventions (Vaida and Blanchard, 2005).

## When the marginal one is available

`loglik_type = "marginal"` returns the value of the outer criterion at
the fit. Its degrees of freedom count the parameters that the marginal
model estimates: the coefficients that are not covered by any penalty,
the dimension of the null space of each penalty without a kink (the
straight line of a smooth), the hyperparameters that the outer search or
a path estimated, and the free parameters of a structural term. A held
hyperparameter is not counted, and a hyperparameter shared by several
terms is counted once. The other coefficients of a penalized term are
integrated out and are not counted.

The value exists only where the outer criterion is
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md) or
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).
Where no outer criterion ran, or where it is
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
there is no marginal likelihood, and requesting one signals an error.

## References

Vaida, F. and Blanchard, S. (2005). Conditional Akaike information for
mixed-effects models. *Biometrika*, 92(2), 351–370.

## See also

[`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
