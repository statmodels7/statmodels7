# The Variance Matrix of a Fit

The variance of the estimated coefficients, over every distribution
parameter's block at once.

## Usage

``` r
# S3 method for class 'StatmodFit'
vcov(
  object,
  type = c("bayesian", "frequentist", "unconditional"),
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

  `"bayesian"`, `"frequentist"` or `"unconditional"`. The first two are
  conditional on the hyperparameters; the third carries their own
  uncertainty as well.

- expected:

  Whether the expected information is used. Defaults to
  [`fit_expected()`](https://statmodels7.github.io/statmodels7/reference/fit_expected.md):
  the expected one where the fit inverted it and the family writes it
  out, the observed Hessian otherwise. A model that carries a structural
  term always uses the observed information.

- approx:

  How the expected information is approximated for a family with no
  closed form: `"opg"` (the default), `"bartlett"`, `"integrate"` or
  `"mc"`. Read only when `expected` is `TRUE` and the family has no
  closed form.

- readable:

  Whether the matrix is over the quantities that the model reports
  (`TRUE`, the default) or over the coordinates on which it was
  estimated (`FALSE`). A structural term declares its quantities, for
  example a persistence, and their variance is obtained by the delta
  method.

- parameter:

  The distribution parameters whose rows and columns are kept. `NULL`
  keeps all of them.

- ...:

  Unused.

## Value

A square matrix with dimnames `parameter:name`. By default it is over
the quantities that the model reports, which include the own parameters
of a structural term. With `readable = FALSE` it is over the
coefficients of every equation, the free parameters of a structural term
following those of its own equation. `parameter` restricts it to the
rows and columns of the named equations. A coordinate that is held or
not identified has `NA` in its row and column.

## Details

**Three matrices.** Writing \\H\\ for the information of the
log-likelihood and \\S\\ for the second derivative of the penalty,
\$\$V_b = (H + S)^{-1}, \qquad V_f = (H+S)^{-1} H (H+S)^{-1}.\$\$ The
first is the posterior variance under the prior whose negative logarithm
is the penalty, and it is the matrix from which an interval around a
penalized term should be built: it carries the smoothing bias as
variance, and for that reason such intervals cover at about their
nominal rate. The second is the sampling variance of the penalized
estimator at a fixed penalty, which is smaller and covers less. With no
penalty \\S = 0\\ and both are \\H^{-1}\\, for a coefficient that the
outer criterion does not estimate. Under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
with the default `marginal = NULL`, in a model without a structural
term, the unpenalized coefficients of every distribution parameter other
than the first (for example the log standard deviation of a Gaussian)
are estimated on the criterion, and `"bayesian"` and `"unconditional"`
take the curvature of that criterion for them. These two types then
differ from `"frequentist"` for those coefficients even without a
penalty.

**A prior that is not Gaussian.** \\V_b\\ is computed in the same way
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

**The first two matrices are conditional on the hyperparameters**, which
they read at the value where the outer search stopped, as though it were
known. That value was estimated from the same data, and the third matrix
adds the variance that this estimation contributes: \$\$V' = V_b + J
V\_\theta J', \qquad J = -(H + S)^{-1} \frac{\partial^2
\rho}{\partial\beta \partial\theta},\$\$ with \\V\_\theta\\ the variance
of the estimated hyperparameters on their own free scale. It is the
delta method applied to the map from the hyperparameter to the penalized
mode (Wood, Pya and Safken, 2016), and the matrix that mgcv returns as
`unconditional = TRUE`. The correction is positive semi-definite, so
\\V'\\ is never narrower than \\V_b\\. It matters where a smoothing
parameter is poorly determined, and it affects the penalized coordinates
much more than the unpenalized ones. It requires the curvature of the
outer criterion, which neither conditional matrix reads, so
`"unconditional"` is more expensive than the other two types.

Where no hyperparameter was estimated by a differentiable criterion
there is nothing to propagate, the correction is exactly zero and
`"unconditional"` returns \\V_b\\ itself. Where a hyperparameter was
estimated and the curvature of the criterion cannot be read at the
fitted point, the conditional matrix is returned with a warning, so that
the narrower matrix is not mistaken for the wider one. Where the
curvature can be read for only some of the hyperparameters, the
correction is a lower bound and a warning states it. See
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md).

**A coefficient that a kinked penalty has set to zero has `NA` in its
row and column.** At zero the penalty is not twice differentiable, so
\\S\\ does not exist there and no curvature can be read. The
coefficients that a lasso or an MCP leaves non-zero have a variance,
which is conditional on that selection;
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
reports this in a note. Away from zero the lasso contributes nothing to
\\S\\, and SCAD and MCP contribute a negative curvature on their concave
part, which makes the variance of a coefficient there larger than that
of the likelihood alone. For these penalties \\V_b\\ is therefore not
the posterior variance under the corresponding prior, because a
curvature read at the mode does not see the mass that the prior puts
near zero.

**Which information.** `expected` selects the matrix used as \\H\\. Its
default is the expected information where the fit inverted it and the
family writes it out in closed form, and the observed Hessian otherwise
([`fit_expected()`](https://statmodels7.github.io/statmodels7/reference/fit_expected.md)).
The two agree asymptotically and not in a sample: one is the information
averaged over the model, the other the curvature of the likelihood at
the data in hand. A model that carries a structural term is always
inverted on the observed information, because neither a filter nor a
mixture over states has an expected information.

Where the family has no closed form, `expected = TRUE` reaches an
approximation, and `approx` selects it. `"opg"`, the default, is the
outer product of the observed scores and costs one gradient;
`"bartlett"` evaluates the expectation itself, a sum over the support
for a discrete family and a quadrature for a continuous one, and is
orders of magnitude more expensive. The expensive route is available and
is not the default.

## References

Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association*, 111(516), 1548–1563.

## See also

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
vapply(c("frequentist", "bayesian", "unconditional"),
       function(ty) sqrt(diag(vcov(fs, type = ty)))[["mu:s(x).z2"]], 0)
#>   frequentist      bayesian unconditional 
#>    0.06507969    0.06551107    0.06623606 
```
