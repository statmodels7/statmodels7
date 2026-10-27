# The Variance Matrix of a Fit

Returns the variance of the estimated coefficients, over every
distribution parameter's block at once.

## Usage

``` r
# S3 method for class 'StatmodFit'
vcov(object, vcov_type = bayesian_vcov(), readable = TRUE, param = NULL, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, which selects the convention of the matrix and, inside the
  object, the information \\H\\.
  [`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
  (the default) returns \\V_b = (H+S)^{-1}\\, the posterior variance
  conditional on the hyperparameters.
  [`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
  returns the sandwich \\V_f = (H+S)^{-1}H(H+S)^{-1}\\, conditional on
  the hyperparameters as well.
  [`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
  returns \\V' = V_b + J V\_\theta J'\\, which adds the uncertainty of
  the hyperparameters. Each takes `information`, which is `NULL` by
  default (the observed information, or the expected information for a
  family with a kink in a parameter, resolved for each fit) and may be
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
  or
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html).
  The conventions are compared at
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md).

- readable:

  Whether the matrix is over the quantities that the model reports
  (`TRUE`, the default) or over the coordinates on which it was
  estimated (`FALSE`). A structural term declares its quantities, for
  example a persistence, and their variance is obtained by the delta
  method.

- param:

  The distribution parameters whose rows and columns are kept, a
  character vector of parameter names. `NULL` keeps all of them.

- ...:

  Unused. Any argument signals an error.

## Value

A square matrix with dimnames `parameter:name`. By default it is over
the quantities that the model reports, which include the own parameters
of a structural term. With `readable = FALSE` it is over the
coefficients of every equation, the free parameters of a structural term
following those of its own equation. `param` restricts it to the rows
and columns of the named equations. A coordinate that is held or not
identified has `NA` in its row and column.

## Details

Writing \\H\\ for the information of the log-likelihood and \\S\\ for
the second derivative of the penalty, \$\$V_b = (H + S)^{-1}, \qquad V_f
= (H+S)^{-1} H (H+S)^{-1}.\$\$ \\V_b\\ is the posterior variance under
the prior whose negative logarithm is the penalty, and it is the matrix
from which an interval around a penalized term should be built: it
carries the smoothing bias as variance, and for that reason such
intervals cover at about their nominal rate. \\V_f\\ is the sampling
variance of the penalized estimator at a fixed penalty, which is smaller
and covers less. With no penalty \\S = 0\\ and both are \\H^{-1}\\, for
a coefficient that the outer criterion does not estimate. Under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
with the default `marginal = NULL`, in a model without a structural
term, the unpenalized coefficients of every distribution parameter other
than the first (for example the log standard deviation of a Gaussian)
are estimated on the criterion, and
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
and
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
take the curvature of that criterion for them. These two conventions
then differ from
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
for those coefficients even without a penalty.

For a prior that is not Gaussian, \\V_b\\ is computed in the same way
for every penalty that is twice differentiable at the mode. \\H + S\\ is
then the curvature of the log-posterior at the mode, and \\V_b\\ is the
Laplace approximation of the posterior variance. The statement about
coverage above holds for a Gaussian prior, whose \\S\\ does not depend
on the coefficients. For another prior \\S\\ is evaluated at the
estimate and can be negative. Under a Student t prior with scale \\s\\
and \\\nu\\ degrees of freedom it is negative for every effect larger
than \\\sqrt{\nu}\\ s\\ in absolute value, so such an effect is shrunk
little and its variance comes mostly from the likelihood. The
approximation is local, and it can be poor where the posterior is skewed
or has more than one mode. Where \\H + S\\ is not positive definite, a
coordinate whose diagonal entry is not positive is held as described
below, and otherwise an error is signaled.

The first two matrices are conditional on the hyperparameters, which are
evaluated at the value where the outer search stopped, as though that
value were known. It was estimated from the same data, and the third
matrix adds the variance that this estimation contributes: \$\$V' =
V_b + J V\_\theta J', \qquad J = -(H + S)^{-1} \frac{\partial^2
\rho}{\partial\beta \partial\theta},\$\$ with \\V\_\theta\\ the variance
of the estimated hyperparameters on their own free scale. It is the
delta method applied to the map from the hyperparameter to the penalized
mode (Wood, Pya and Safken, 2016), and the matrix that mgcv returns as
`unconditional = TRUE`. The correction is positive semi-definite, so
\\V'\\ is never narrower than \\V_b\\. It matters where a smoothing
parameter is poorly determined, and it affects the penalized coordinates
much more than the unpenalized ones. It requires the curvature of the
outer criterion, which the two conditional matrices do not read, so
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
is more expensive than the other two conventions.

Where no hyperparameter was estimated by a differentiable criterion
there is nothing to propagate, the correction is exactly zero and
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
returns \\V_b\\ itself. Where a hyperparameter was estimated and the
curvature of the criterion cannot be read at the fitted point, the
conditional matrix is returned with a warning, so that the narrower
matrix is not mistaken for the wider one. Where the curvature can be
read for only some of the hyperparameters, the correction is a lower
bound and a warning states it. See
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md).

A coefficient that a kinked penalty has set to zero has `NA` in its row
and column. At zero the penalty is not twice differentiable, so \\S\\
does not exist there and the curvature cannot be read. The coefficients
that a lasso or an MCP leaves non-zero have a variance, which is
conditional on that selection;
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
reports this in a note. Away from zero the lasso contributes nothing to
\\S\\, and SCAD and MCP contribute a negative curvature on their concave
part, which makes the variance of a coefficient there larger than that
of the likelihood alone. For these penalties \\V_b\\ is therefore not
the posterior variance under the corresponding prior, because a
curvature read at the mode does not include the mass that the prior puts
near zero.

The coefficients that the marginal criterion of
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimates take the curvature of that criterion in the Bayesian and the
unconditional matrices, and the criterion reads the information of its
own `information` argument. The `information` of `vcov_type` selects the
matrix used as \\H\\ everywhere else, in the sandwich of the frequentist
matrix and in the penalized block of the other two.

The `information` of `vcov_type` selects the matrix used as \\H\\. Its
default is `NULL`, which is resolved for each fit by the family of the
fit with the rule of
[`distributions7::vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.html):
the matrix is
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html),
the curvature of the likelihood at the data in hand, when the
log-likelihood is twice differentiable in every parameter of the family
(every
[`distributions7::params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.html)
is at least 2), and
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
with `method = NULL`, for the whole matrix, otherwise. The choice does
not depend on the information with which the fit was computed
(`inner_information` of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)).
An object in `information` overrides the rule.
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
selects the Fisher information, the information averaged over the model.
The two agree asymptotically and not in a sample. A model that carries a
structural term is always inverted on the observed information, because
a filter and a mixture over states do not have an expected information.
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
states the matrix when it is not the observed one.

The observed information is singular where the second derivative of the
log-likelihood in a coefficient is zero at the estimates, as for the
location of a Laplace response. The family has a kink in that parameter,
so the default `NULL` reads the expected information, which is finite
there, and the variance is available. With
`information = observed_info()` the rows and columns of the coordinates
that it does not identify are missing, as described above for a flat
direction.

The `method` of
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
selects how the expected information is obtained. With `method = NULL`
the closed form of the family is used where the family has one, and
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
the outer product of the observed scores, which costs one gradient,
otherwise.
[`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
evaluates the expectation itself, a sum over the support for a discrete
family and a quadrature for a continuous one, and is orders of magnitude
more expensive;
[`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html)
simulates it.
[distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html)
compares the three methods.

## References

Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association*, 111(516), 1548–1563.

## See also

[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the three conventions side by side,
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md),
which builds the third matrix's correction

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(80))
dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
sqrt(diag(vcov(fit)))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>        0.08770104        0.14865515        0.08006408 

# With a penalized term the three differ for a penalized coefficient, and
# the widest is the one that does not read the smoothing parameter as known.
ds <- data.frame(x = runif(200))
ds$y <- sin(2 * pi * ds$x) + rnorm(200, sd = 0.3)
fs <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), ds)
vapply(list(frequentist_vcov(), bayesian_vcov(), unconditional_vcov()),
       function(ty) sqrt(diag(vcov(fs, vcov_type = ty)))[["mu:s(x).z2"]], 0)
#> [1] 0.06510115 0.06553290 0.06625519

# the information is chosen inside the convention
sqrt(diag(vcov(fit, vcov_type = bayesian_vcov(
  information = distributions7::expected_info()))))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>        0.08770104        0.14865515        0.08006408 
```
