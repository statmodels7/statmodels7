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

**Three matrices, and they differ only when something is penalized.**
Writing \\H\\ for the information of the log-likelihood and \\S\\ for
the second derivative of the penalty, \$\$V_b = (H + S)^{-1}, \qquad V_f
= (H+S)^{-1} H (H+S)^{-1}.\$\$ The first is the posterior variance under
the prior the penalty is the negative logarithm of, and it is what an
interval around a penalized term should be built from: it carries the
smoothing bias as though it were variance, and that is why such
intervals cover at about their nominal rate. The second is the sampling
variance of the penalized estimator at a fixed penalty, which is smaller
and covers less. With no penalty \\S = 0\\ and both are \\H^{-1}\\.

**Both of those are conditional on the hyperparameters**, read at the
value the outer search stopped at as though it had been known. It was
estimated from the same data, and the third matrix adds what that costs:
\$\$V' = V_b + J V\_\theta J', \qquad J = -(H + S)^{-1} \frac{\partial^2
\rho}{\partial\beta \partial\theta},\$\$ with \\V\_\theta\\ the variance
of the estimated hyperparameters on their own free scale. It is the
delta method applied to the map from the hyperparameter to the penalized
mode (Wood, Pya and Safken, 2016), and the matrix mgcv returns as
`unconditional = TRUE`. The correction is positive semi-definite, so
\\V'\\ is never narrower than \\V_b\\; measured on a univariate smooth
it widens the band of the fitted mean by 1.1 per cent on average and 7.6
per cent at its widest point at \\n = 200\\, and by 0.2 and 1.3 per cent
at \\n = 2000\\. It is where a smoothing parameter is poorly determined
that it matters: the coordinates a penalty compresses widen by as much
as 88 per cent on a fit whose \\\log\lambda\\ carries a standard
deviation of 2.4, while the unpenalized coordinates beside them move in
the sixth decimal.

It costs four to five times the matrix it is added to – 2.7 ms against
10.6 ms at \\n = 200\\ and 4.0 against 20.0 at \\n = 2000\\ – because it
reads the outer criterion's own curvature, which neither conditional
matrix asks for.

Where no hyperparameter was estimated by a differentiable criterion
there is nothing to propagate, the correction is exactly zero and
`"unconditional"` returns \\V_b\\ itself. Where one was estimated and
its curvature cannot be read – a shared hyperparameter, or one the
search left at the edge of its range – the conditional matrix is
returned WITH A WARNING, since a reader who asked for the wider matrix
and silently received the narrower one would report the wrong thing. See
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md).

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
#>        0.08659784        0.14678520        0.07905694 

# With a penalized term the three differ, and the widest is the one that
# does not read the smoothing parameter as known.
ds <- data.frame(x = runif(200))
ds$y <- sin(2 * pi * ds$x) + rnorm(200, sd = 0.3)
fs <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds)
vapply(c("frequentist", "bayesian", "unconditional"),
       function(ty) sqrt(diag(vcov(fs, type = ty)))[[1L]], 0)
#>   frequentist      bayesian unconditional 
#>    0.02085794    0.02085794    0.02085794 
```
