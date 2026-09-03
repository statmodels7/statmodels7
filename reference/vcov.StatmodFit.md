# The Variance Matrix of a Fit

The variance of the estimated coefficients, over every distribution
parameter's block at once.

## Usage

``` r
# S3 method for class 'StatmodFit'
vcov(
  object,
  type = c("bayesian", "frequentist"),
  expected = NULL,
  approx = c("opg", "bartlett", "integrate", "mc"),
  readable = TRUE,
  parameter = NULL,
  ...
)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- type:

  `"bayesian"` or `"frequentist"`.

- expected:

  Whether the expected information is used. Defaults to
  [`fit_expected()`](https://statmodels7.github.io/statmodels7/reference/fit_expected.md):
  the expected one where the fit inverted it and the family writes it
  out, the observed Hessian otherwise.

- approx:

  How the expected information is approximated for a family with no
  closed form: `"opg"` (the default), `"bartlett"`, `"integrate"` or
  `"mc"`. Read only when `expected` is `TRUE` and the family has no
  closed form.

- ...:

  Unused.

## Value

A square matrix over the stacked coefficients, with dimnames
`parameter:coefficient`.

## Details

**Two matrices, and they differ only when something is penalized.**
Writing \\H\\ for the information of the log-likelihood and \\S\\ for
the second derivative of the penalty, \$\$V_b = (H + S)^{-1}, \qquad V_f
= (H+S)^{-1} H (H+S)^{-1}.\$\$ The first is the posterior variance under
the prior the penalty is the negative logarithm of, and it is what an
interval around a penalized term should be built from: it carries the
smoothing bias as though it were variance, and that is why such
intervals cover at about their nominal rate. The second is the sampling
variance of the penalized estimator at a fixed penalty, which is smaller
and covers less. With no penalty \\S = 0\\ and both are \\H^{-1}\\.

**A coefficient a kinked penalty has set to zero has no row.** At zero
the penalty is not twice differentiable, so \\S\\ does not exist there
and no curvature can be read; the entry is `NA`. The coefficients a
lasso or an MCP left non-zero do get a variance, and it is conditional
on that selection, which
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
says in a note instead of leaving a reader to assume otherwise. **Which
information.** `expected` says which matrix \\H\\ is. Its default is the
expected information where the fit inverted it AND the family writes it
out in closed form, and the observed Hessian otherwise
([`fit_expected()`](https://statmodels7.github.io/statmodels7/reference/fit_expected.md)).
The two agree asymptotically and not in a sample: one is the information
averaged over the model, the other the curvature of the likelihood at
the data in hand.

Where the family has no closed form, `expected = TRUE` reaches an
approximation and `approx` says which. `"opg"`, the default, is the
outer product of the observed scores and costs one gradient;
`"bartlett"` evaluates the expectation itself, a sum over the support
for a discrete family and a quadrature for a continuous one, and is
orders of magnitude dearer – 89.06 s against 0.64 s on a Poisson-inverse
gaussian regression at \\n = 500\\. The expensive route is reachable and
is not the default.

## See also

[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(80))
dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
sqrt(diag(vcov(fit)))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>        0.08659784        0.14678520        0.07905694 
```
