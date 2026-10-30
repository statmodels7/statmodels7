# Summarize a Fitted Model

Returns one coefficient table per distribution parameter, carrying the
estimate, its standard error, the statistic of the test named by `test`
(Wald's by default), its p-value and an interval, together with the
degrees of freedom, the information criteria and the notes that qualify
the numbers.

## Usage

``` r
# S3 method for class 'StatmodFit'
summary(
  object,
  level = 0.95,
  vcov_type = bayesian_vcov(),
  correct = FALSE,
  test = wald_test(),
  max_coef = NULL,
  ...
)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- level:

  The confidence level.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, the variance matrix, passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md):
  [`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
  (the default),
  [`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
  or
  [`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md).
  The information from which the standard errors are computed is an
  argument of each of them
  ([`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
  by default). The last carries the uncertainty of the hyperparameters
  into every standard error under it, where `correct` carries the same
  uncertainty into the degrees of freedom; the two are the same quantity
  read on two surfaces, and both can be used in one call. The
  conventions are compared at
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md).
  The effective degrees of freedom of the summary are read on this same
  information (see the section on the degrees of freedom).

- correct:

  Whether the degrees of freedom include the cost of estimating the
  hyperparameters. The ordinary count treats them as known although they
  were estimated from the same data, so a criterion built on it is too
  generous. See
  [`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md).
  The default is `FALSE`, which leaves the count unchanged; the
  correction is zero where no hyperparameter was estimated by a marginal
  criterion.

- test:

  A
  [test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
  object, the statistic that the coefficient tables report.
  [`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md)
  (the default) reads the variance matrix and refits nothing.
  [`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
  [`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
  and
  [`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md)
  cost one restricted refit per row and are reported as a signed root.
  The four are compared at
  [test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md);
  see also
  [`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md).
  A character string is rejected with an error.

- max_coef:

  How many coefficient rows each block of the printed summary shows,
  `Inf` or `NA` for all of them. `NULL`, the default, reads the option
  `statmodels7.summary_max_coef`, and 10 where that is unset. A
  hyperparameter row is never hidden, whatever this is, and a block of
  twelve rows or fewer is never abridged. The value is carried on the
  result, so `summary(fit, max_coef = Inf)` prints in full at the
  console;
  [`print.StatmodSummary()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodSummary.md)
  takes the same argument for an object already in hand.

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

A
[`StatmodSummary()`](https://statmodels7.github.io/statmodels7/reference/StatmodSummary-class.md).

## Details

Each distribution parameter is read as blocks, not as one list of
coefficients, because most coefficients of a fitted model are not
quantities of interest in themselves. The blocks are

- the parametric terms:

  every unpenalized term together, one row per coefficient, which is the
  ordinary table.

- one block per smooth:

  the linear component's coefficient where the construction carries one,
  the smoothing parameter, and the effective degrees of freedom. The
  coefficients of the wiggly part are not shown: they are coordinates in
  an orthonormal basis, a single one of them has no interpretation, and
  together they determine the edf.

- one block per random effect:

  the parameters of the distribution of the effects, usually called the
  variance components, the edf and a summary of the predicted effects.
  The effects themselves, of which there is one per level, are not
  listed.

- one block per selection:

  a lasso, a SCAD or an MCP: its hyperparameters, how many coefficients
  survived, and those coefficients. The ones set exactly to zero are
  counted, not listed.

- one block per other penalized term:

  its coefficients, which stay interpretable under a ridge, together
  with its hyperparameters.

- one block per break-point term:

  the quantities that the term reports (slopes, changes and positions),
  with standard errors from the delta method applied to the coefficients
  through which the term is fitted. A position carries no test.

- one block per structural term:

  the term's own parameters (a filter's level, loadings and
  persistences, or the transition probabilities of a regime chain), with
  a level that an intercept holds marked as held.

## The effective degrees of freedom

The edf of a term is its share of \\\mathrm{tr}\[(H+S)^{-1}H\]\\, with
\\H\\ the information of the log-likelihood and \\S\\ the second
derivative of the penalty, read over the coefficients of every equation
at once. The summary reads \\H\\ on the information of `vcov_type`
([`distributions7::default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.html)
of the family where the `information` of the variance object is `NULL`),
which is the matrix that the standard errors above the edf invert. A
call with `vcov_type = bayesian_vcov(information = expected_info())`
therefore prints the standard errors and the edf on the expected
information. The count stored in `fit@edf` was read at fit time on the
information of the outer criterion: the observed information by default
under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).
With the defaults the two coincide, and the summary recomputes the count
only where `vcov_type` names the other matrix.

Three cases keep the count stored in the fit. Under
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
the count is part of the criterion that prices the fit, and it is read
on the `information` of that criterion, so the summary leaves it as it
is and the printed criterion stays reproducible from it. A model with a
structural term is read on the joint vector of the coefficients and the
term's own parameters, as before. A variance object without an
`information` property has no matrix to read.

The eigenvalues of \\(H+S)^{-1}H\\ are \\v^\top H v/(v^\top H v + v^\top
S v)\\, for the generalized eigenvectors \\v\\. They lie in \\\[0,
1\]\\, and the total in \\\[0, p\]\\ with \\p\\ the number of
coefficients, when \\H\\ is positive semidefinite. The observed
information can be indefinite (a Student t at an outlier, for example).
A direction with \\v^\top H v \< 0\\ then contributes a negative share,
an edf can be underestimated, and the bounds above can fail. The summary
adds a note when the observed information used for the count has an
eigenvalue below \\-p\\\epsilon\\\lambda\_{\max}\\ (\\\epsilon\\ the
machine epsilon and \\\lambda\_{\max}\\ the largest eigenvalue), the
floor with which the package separates a flat direction from rounding.
The expected information is positive semidefinite and carries no such
note.

A hyperparameter is the first row of its block, since it governs every
coefficient under it, and its name carries a mark in square brackets
that gives the source of the value. A hyperparameter estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
maximizes a twice differentiable criterion, so it carries a standard
error and an interval, both read on the free scale that its link defines
and mapped back
([`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)).
One chosen by
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
along a path over a kinked penalty is the argument of a minimum over a
grid, so the row names the criterion and leaves the remaining columns
empty: there is no curvature at such a point from which to read a
standard error. A smooth penalty under
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) is
minimized over its free scale, and its row also has no standard error,
because the criterion is not a likelihood. A hyperparameter that the
term holds is marked fixed, and one that stayed at the starting value of
its term because the outer criterion was not run is marked start.

The statistic column reports the test named by `test`. By default it is
Wald's, the estimate over its standard error, which is the only one of
the four that needs no refit and is read only from the unrestricted fit.
The other three – the likelihood ratio, Rao's score and Terrell's
gradient – hold the coefficient at the null value and refit, one refit
per row, and are reported as the signed root of their \\\chi^2_1\\
statistic with the sign of \\\hat\beta_j - b\\, so that the column
carries one kind of quantity whichever test was requested and Wald's own
\\z\\ is the case where that root is exact. The p-value is the
\\\chi^2_1\\ one in every case, and for Wald the two readings agree
identically.

Only the statistic and its p-value change with `test`: the interval
stays the Wald one, since inverting a test requires many refits for each
coefficient.
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
inverts the same test through its `test` argument, one coefficient at a
time. A row for which the restricted fit is not defined – a coefficient
under a kinked penalty, an aliased one, any coefficient of a model
carrying a structural term – reports `NA` instead of the Wald value, so
the column never carries two tests at once. Every other row is tested,
including a smooth's own coordinates and a random effect's, on the
penalized objective and with the reading described in
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md).
A quantity of a break-point term that is one coefficient carries the
restricted statistic of that coefficient; a quantity that combines
several coefficients, and the position of a break-point, carry none.

The meaning of a Wald p-value depends on the row. For an unpenalized
coefficient it is the ordinary p-value. For a coefficient in a penalized
block it is conditional on the smoothing parameter, which was not
estimated jointly with it, and it does not account for the shrinkage of
the estimate toward zero. For a coefficient that a kinked penalty
selected, the row exists only because that coefficient survived the
selection, and a naive interval there under-covers.

The degrees of freedom are the effective ones, summed over the terms, so
a penalized term counts the degrees of freedom that it uses instead of
the number of its columns. The information criteria are built on that
count: `aic` is \\-2\ell + 2\\\mathrm{df}\\ and `bic` replaces the
factor 2 by the logarithm of the sum of the prior weights, which are
frequency weights.

## Reading the foot

Below the criteria the summary describes the point at which the fit
stopped with a few numbers from
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md).
The summary prints no verdict. Each number is dimensionless, so it does
not change with the units of a coefficient.

- `inner max |grad|/se`:

  the largest \\\|g_j\|/\sqrt{A\_{jj}}\\ over the coefficients, where
  \\g\\ is the gradient of the penalized objective and \\A\\ its
  curvature, the penalized information. \\1/\sqrt{A\_{jj}}\\ is the
  standard error of coefficient \\j\\ conditional on the others, which
  is smaller than the marginal one that the tables print, and the ratio
  equals \\\sqrt{2\Delta_j}\\ with \\\Delta_j = g_j^2/(2A\_{jj})\\ the
  gain in log-likelihood of a Newton step in that coordinate alone. At a
  mode it is small. For a model carrying a score-driven filter it covers
  the filter's own parameters as well. It leaves out a coefficient that
  a kinked penalty has set to zero, a frozen working block, and an
  aliased column, and it is read on the directions that the kinks of a
  break-point sitting on an observation leave free.

- `inner min eigen`:

  the smallest eigenvalue of \\D^{-1/2}AD^{-1/2}\\ with \\D =
  \mathrm{diag}(A)\\, the curvature rescaled to a unit diagonal. It is
  at most one. Near one the coordinates are well separated; near zero
  some combination of them is barely identified, as two nearly collinear
  covariates or an intercept beside a random effect make it; below zero
  the point is not a mode.

- `zeros max |score|/kink`:

  printed only where a kinked penalty (a lasso, an elastic net, a SCAD
  or an MCP) set coefficients exactly to zero, together with how many.
  At a zero the penalty has no derivative and the coefficient no
  standard error, so the condition for an optimum is \\\|s_j\| \le
  \kappa_j\\: the derivative of the log-likelihood in that coefficient
  within the size of the kink, \\\lambda\\ for a lasso. The number is
  the largest \\\|s_j\|/\kappa_j\\. At or below one every zero is an
  optimum; above one a coefficient would improve the objective by
  leaving zero. A reading \\r\\ below one means that the coefficient
  with the largest ratio is held by the kink with a margin of \\1 - r\\;
  if the other coefficients were fixed it would leave zero at a
  \\\lambda\\ smaller by the factor \\r\\. Along a path the others move
  too, and the threshold moves with them. The number is not a test of
  significance, only a check on whether the zeros agree with the
  \\\lambda\\ of the fit, and values close to one are ordinary where
  \\\lambda\\ was chosen on a grid, because the next coefficient to
  enter is close to its threshold by construction.

- `outer max |grad|/se` and `min eigen`:

  the same two readings for the criterion that estimated the
  hyperparameters, on the free scale that their links define, printed
  only where
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  or
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  ran. The curvature is the negated Hessian of the criterion. A
  hyperparameter that ran to the edge of its chart is left out, because
  its curvature collapses with its gradient; where every one did the
  line reads `not available: at a boundary`. A hyperparameter chosen
  along a path has no gradient, and there is no outer line.

## See also

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md),
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md),
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md),
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md),
which computes the same four statistics for one coefficient against a
value that need not be zero, and
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(120))
dd$y <- 1 + 2 * dd$x + rnorm(120, sd = 0.4)
summary(statmod(y ~ x | sigma ~ x,
                distributions7::gaussian1_distrib(), dd))
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ x | sigma ~ x, distrib = distributions7::gaussian1_distrib(), 
#>             data = dd)
#> 
#> Distribution: gaussian1     Observations: 120
#> 
#> === mu   [identity link]
#> 
#> Parametric terms
#>                estimate      se     z       p  lower upper
#>   (Intercept)     1.063 0.07583 14.02 < 1e-16 0.9143 1.212
#>   x               1.903  0.1313 14.49 < 1e-16  1.645  2.16
#> 
#> === sigma   [log link]
#> 
#> Parametric terms
#>                estimate     se       z         p   lower   upper
#>   (Intercept)   -0.9559 0.1338  -7.144 9.051e-13  -1.218 -0.6936
#>   x             0.02082 0.2277 0.09145    0.9271 -0.4254   0.467
#> 
#> 95% intervals, bayesian variance
#> conditional log-likelihood -55.852571    effective df 4.00
#> cAIC 119.705    cBIC 130.855
#> fitted in 114 ms
#> inner   max |grad|/se 1.6e-15   min eigen 0.12
#> outer   max |grad|/se 1.8e-14   min eigen 0.13
```
