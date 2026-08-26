# The Maximized Log-Likelihood of a Fit

The value R's convention expects, carrying the degrees of freedom and
the number of observations.
[`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
is the other thing: the model evaluated at parameters and data of the
caller's choosing.

## Usage

``` r
# S3 method for class 'StatmodFit'
logLik(object, type = c("conditional", "marginal"), ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- type:

  `"conditional"` (the default) or `"marginal"`.

- ...:

  Unused.

## Value

A `logLik` object: a single number with attributes `df`, the degrees of
freedom, and `nobs`. [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html)
and [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) read it, so
`AIC(fit)` is the conditional AIC.

## Which likelihood, and why it matters

The default is the **conditional** one: the log-density at the fitted
coefficients, a penalized coefficient among them, paired with the
effective degrees of freedom \\\mathrm{tr}\[(H+S)^{-1}H\]\\. A criterion
built on that pair is the conditional AIC, and it asks how well the
model describes the groups, curves and states actually observed.

A mixed-model package reports the **marginal** likelihood instead: the
random effects are integrated out and the count is the number of
estimated parameters, variance components among them. That asks about
the population, and its AIC is a different number, not comparable with
the conditional one.

Neither convention allows the halves to be mixed: a marginal likelihood
against an effective count, or a conditional one against a parameter
count (Vaida and Blanchard, 2005).

Measured on a random intercept over 120 groups, the two readings of the
same fit are -4148.59 on 115.72 effective degrees of freedom and
-4372.79 on 4, against `lme4::lmer`'s marginal -4371.71 on 4.

## When the marginal one is available

`type = "marginal"` returns the value the outer criterion evaluated
while choosing the hyperparameters, with the number of estimated
parameters as its degrees of freedom. It exists only where a marginal
criterion actually ran, which is
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md) or
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).

Where the hyperparameters were held, or chosen by
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
there is no marginal likelihood to report and asking for one signals an
error instead of returning a number that would look like one.

## References

Vaida, F. and Blanchard, S. (2005). Conditional Akaike information for
mixed-effects models. *Biometrika*, 92(2), 351–370.

## See also

[`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
