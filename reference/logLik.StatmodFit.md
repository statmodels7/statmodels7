# The Maximized Log-Likelihood of a Fit

The value R's convention expects, carrying the degrees of freedom and
the number of observations.
[`loglik`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
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
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- type:

  `"conditional"` (default) or `"marginal"`.

- ...:

  Unused.

## Value

A `logLik` object.

## Details

**Which likelihood, and it matters.** The default is the CONDITIONAL
one: the log-density at the fitted coefficients, a penalized coefficient
among them, paired with the effective degrees of freedom
\\\mathrm{tr}\[(H+S)^{-1}H\]\\. A criterion built on the pair is the
conditional AIC, and its question is how well the model describes the
groups, curves and states actually observed.

A mixed-model package reports the MARGINAL likelihood instead: the
random effects are integrated out and the count is the number of
estimated parameters, variance components among them. Its question is
about the population, and its AIC is a different number that is not
comparable with the conditional one. Mixing the two – a marginal
likelihood against an effective count, or the reverse – is what neither
convention allows (Vaida and Blanchard, 2005).

`type = "marginal"` returns the value the outer criterion evaluated
while choosing the hyperparameters, with the number of estimated
parameters as its degrees of freedom. It is available only where a
marginal criterion actually ran:
[`ml`](https://statmodels7.github.io/statmodels7/reference/reml.md) or
[`reml`](https://statmodels7.github.io/statmodels7/reference/reml.md).
Where the hyperparameters were held, or found by a prediction criterion,
there is no marginal likelihood to report and asking for one is an error
rather than a number that would look like one.

## References

Vaida, F. and Blanchard, S. (2005). Conditional Akaike information for
mixed-effects models. *Biometrika*, 92(2), 351–370.

## See also

[`loglik`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
