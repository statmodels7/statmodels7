# Information Matrices Available to a Fit

An information object selects the matrix of second derivatives of the
log-likelihood that a computation inverts or takes the determinant of.
Three objects are available:
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
and
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html),
which distributions7 defines and documents, and
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md),
which statmodels7 defines. All three inherit from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html).

Three arguments take them, and each one selects the matrix for its own
computation, so that the choices are independent:

- `inner_information` of
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md):

  The matrix of the inner fit, which solves for the coefficients at
  fixed hyperparameters. With
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  the matrix is inverted at each step; with an optimizers7 optimizer it
  is the Hessian that the optimizer receives. The default is
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md).

- `information` of
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md):

  The matrix of the outer criterion: its determinant enters the Laplace
  approximation of
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  and
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  and its trace gives the effective degrees of freedom of
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md).
  The effective degrees of freedom stored in the fitted model
  (`fit@edf`) are read on this matrix whenever an outer criterion ran,
  and on the matrix of `inner_information` only for a model without one.
  The default is
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html).

- `information` of
  [`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md),
  [`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
  and
  [`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md):

  The matrix inverted for the standard errors. The default is `NULL`,
  which selects the matrix by the family of the fit: the observed
  information where the log-likelihood is twice differentiable in every
  parameter of the family, and
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  otherwise (see "The default of the variance objects" below).
  [`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
  reads the effective degrees of freedom it prints on this matrix too,
  so that they are shares of the same \\H\\ as the standard errors
  beside them.

## Value

No value. The topic describes the objects returned by
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html),
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
and
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md).

## Details

Let \\\ell(\beta)\\ be the log-likelihood of the coefficients \\\beta\\.

- [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html):

  The negative Hessian of the log-likelihood at the data, \\-\partial^2
  \ell / \partial\beta\\\partial\beta^\top\\. It is the curvature of the
  likelihood in the sample at hand, and every family has it. It may be
  indefinite far from the mode, and for a family such as the Student t
  at an outlier. It is always the matrix of the data, and the expected
  information is never substituted for it. A step that needs a positive
  definite matrix floors the eigenvalues of the penalized matrix (see
  [`solve_chol()`](https://statmodels7.github.io/statmodels7/reference/solve_chol.md)).
  The observed information of a coefficient whose second derivative is
  zero at the estimates, as for the location of a Laplace response, is
  singular, and a variance that inverts it is not available there.

- [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html):

  The Fisher information, \\\mathbb{E}\[-\partial^2 \ell /
  \partial\beta\\\partial\beta^\top\]\\, a function of the parameters
  alone. Its `method` selects how the expectation is obtained (see "The
  method of the expected information" below).

- [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md):

  A choice that is settled once the family is known. The rule depends on
  where the object is passed and is given below. The fitted model
  records the settled choice.

## Which function accepts which object

|  |  |  |  |  |
|----|----|----|----|----|
| Argument | observed | expected | auto | `NULL` |
| `inner_information` of [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md) | yes | yes | yes (default) | no |
| `information` of [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md) and [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md) | yes (default) | yes | yes | no |
| `information` of [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) and [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) | yes (default) | yes | no | no |
| `information` of the three variance objects | yes | yes | no | yes (default) |

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and the variance objects reject
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
because the quantity that they compute is defined on one named matrix.

## The method of the expected information

The `method` of
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
and of
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
selects how the expected information, and the derivatives of it that an
exact outer gradient reads, are obtained. It is `NULL` or one of
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
[`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
and
[`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html);
the three are compared in
[distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html).

- With `method = NULL`, the default, the closed form of the family is
  used where the family writes the expected information out (see
  [`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)),
  and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise. The values are those of the closed form for every shipped
  family, and the outer criterion under the default is the same
  criterion as before the argument existed.

- An object is always used, also for a family that has the closed form.

- The three arguments above are independent, and so are their methods.

In an outer criterion
([`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md))
on the expected information, one method object defines the matrix and
every derivative of it that the exact outer gradient and Hessian read.
The derivative is the derivative of the matrix that the method defines:

- for the closed form, the closed form of the derivative;

- for
  [`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
  and
  [`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html),
  the derivative of the expectation, which moves the density as well as
  the integrand;

- for
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html),
  the exact derivative of the outer product of the scores at the data,
  which treats the data as fixed.

Where the method is `NULL` and the family has no closed form for the
derivative, the matrix and its derivatives are both
[`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html).
The outer product of the scores depends on the response, so a criterion
built on it is a smooth function of the hyperparameters with exact
derivatives, and it is not the Laplace approximation of an expectation:
the matrix is an estimate of the expected information, noisier than a
computed expectation and about 25 times cheaper. Measured on a
density-only family (n = 300, ten basis functions),
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
took 0.8 s with `method = NULL` and 25.3 s with
`method = expect_integrate(identity = "score")`, and the smoothing
parameter was 12.57 and 13.00, against 13.02 on the observed
information. The computed expectation is selected with
`method = expect_integrate(identity = "score")`.
[`distributions7::expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.html)
draws a new sample at each call of the matrix and of each derivative, so
that the derivative is that of the expectation and not of the sampled
matrix; the outer search is then noisy at the level of the Monte Carlo
error.

## The rule of auto_info()

In `inner_information`, the object is settled by the family and, with
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
by the run.

- The expected information is chosen where the family writes it exactly
  (see
  [`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html))
  and at an acceptable cost (see
  [`distributions7::expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.html)).
  The observed information is chosen for every other family.

- With
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  on the expected information, the run continues on the observed
  information as soon as Fisher scoring contracts the score slowly (see
  `switch_after` and `switch_rate` in
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)).

- With
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  on the observed information, the expected information takes the step
  where the observed step, after the repair of an indefinite matrix,
  finds no acceptable point. The `repair_floor` of
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  sets the relative eigenvalue under which the matrix counts as
  indefinite.

- With an optimizers7 optimizer, the settled matrix is the Hessian that
  the optimizer receives (the argument `he` of
  [`optimizers7::minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.html)).
  Neither the switch nor the fallback exists there, because both belong
  to the loop of
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md).

In `information` of
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
the expected information is chosen for a Student t, a skew t or a Cauchy
response (wrapped or not), whose observed information is not positive
definite at an outlier, and the observed information otherwise.

## The default of the variance objects

The `information` of
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
and
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
is `NULL`. It is resolved for each fit by the family of the fit, in the
same way as the `information` of
[`distributions7::vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.html):
the matrix is
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
when
[`distributions7::params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.html)
is at least 2 for every parameter of the family, that is, when the
log-likelihood is twice differentiable in every parameter, and
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
with `method = NULL` otherwise, for the whole matrix. A family with a
kink in a parameter (the Laplace, the elastic net) has a singular or
undefined observed information in that parameter, and the closed-form
expected information gives finite standard errors.
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
states the matrix whenever it is not the observed information, and a
printed variance object states the information that it names.

A model that carries a structural term (a score-driven filter or a
mixture over states) is fitted and inverted on the observed information
whatever the object says, because those terms have no expected
information.

## See also

[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
and
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
for the two matrices,
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
for the third object,
[AutoInfo](https://statmodels7.github.io/statmodels7/reference/AutoInfo.md)
for its class,
[distributions7::expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.html)
for the methods of the expected information,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the functions that take an object, and
[iwls_solvers](https://statmodels7.github.io/statmodels7/reference/iwls_solvers.md)
for the factorization of the matrix.

## Examples

``` r
distributions7::observed_info()
#> Observed information
distributions7::expected_info()
#> Expected information (closed form where the family has one, otherwise the outer product of the scores)
auto_info(method = distributions7::expect_integrate(identity = "score"))
#> Information chosen against the family (expectation by integration against the density (identity = score))

set.seed(1)
dd <- data.frame(x = runif(120, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
G <- distributions7::gaussian1_distrib()
fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
               inner_information = distributions7::expected_info(),
               outer_criterion = reml(
                 information = distributions7::observed_info()))
fit@methods$smooth
#> iwls: qr, expected information
#>   maxit 1000, tol 1e-06

# The standard errors are chosen after the fit, independently of both.
sqrt(diag(vcov(fit, vcov_type = bayesian_vcov(
  information = distributions7::expected_info()))))
#>    mu:(Intercept)       mu:s(x).lin        mu:s(x).z1        mu:s(x).z2 
#>        0.02545228        0.02555900        0.02780178        0.07264020 
#>        mu:s(x).z3        mu:s(x).z4        mu:s(x).z5        mu:s(x).z6 
#>        0.14724547        0.22451921        0.48824526        0.51072358 
#> sigma:(Intercept) 
#>        0.06684662 
```
