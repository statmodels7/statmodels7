# Summarize a Fitted Model

One coefficient table per distribution parameter, carrying the estimate,
its standard error, the Wald statistic, its p-value and an interval,
together with the degrees of freedom, the information criteria and the
qualifications the numbers carry.

## Usage

``` r
# S3 method for class 'StatmodFit'
summary(
  object,
  level = 0.95,
  type = c("bayesian", "frequentist", "unconditional"),
  expected = NULL,
  approx = c("opg", "bartlett", "integrate", "mc"),
  correct = FALSE,
  test = c("wald", "lr", "score", "gradient"),
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

- type:

  Which variance matrix: passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).
  The `"unconditional"` one carries the hyperparameters' own uncertainty
  into every standard error under it, where `correct` carries the same
  uncertainty into the degrees of freedom; the two are the same quantity
  read on two surfaces and can be asked for together.

- expected:

  Which information the standard errors are read off, passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).
  `NULL`, the default, takes the expected one where the fit inverted it
  and the family writes it out, and the observed Hessian otherwise;
  `TRUE` or `FALSE` names one.

- approx:

  How the expected information is approximated for a family with no
  closed form, passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md):
  `"opg"` (the default), `"bartlett"`, `"integrate"` or `"mc"`.

- correct:

  Whether the degrees of freedom carry what the estimation of the
  hyperparameters cost. The ordinary count reads them as known, and they
  were chosen from the same data, so a criterion built on it is too
  generous. See
  [`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md).
  Defaults to `FALSE` because it changes a number a reader may be
  comparing with an earlier fit; it is zero where no hyperparameter was
  estimated.

- test:

  Which statistic the coefficient tables report: `"wald"`, the default,
  or `"lr"`, `"score"` or `"gradient"`, each of which costs one
  restricted refit per row. See
  [`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md).

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

**Each distribution parameter is read as blocks, not as one list of
coefficients**, because most of a fitted model's coefficients are not
quantities anybody reads. The blocks are

- the parametric terms:

  every unpenalized term together, one row per coefficient, which is the
  ordinary table.

- one block per smooth:

  the linear component's coefficient where the construction carries one,
  the smoothing parameter, and the effective degrees of freedom. The
  coefficients of the wiggly part are not shown: they are coordinates in
  an orthonormal basis and say nothing one at a time, while what they
  say together is the edf.

- one block per random effect:

  the parameters of the effects' distribution, usually called the
  variance component, and the edf. Not the effects themselves, of which
  there is one per level.

- one block per selection:

  a lasso, a SCAD or an MCP: its hyperparameters, how many coefficients
  survived, and those coefficients. The ones set exactly to zero are
  counted, not listed.

- one block per other penalized term:

  its coefficients, which stay interpretable under a ridge, together
  with its hyperparameters.

**A hyperparameter is the first row of its block**, since it governs
every coefficient under it, and the cell where its standard error would
be says what put the value there. One estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
maximizes a twice differentiable criterion, so it carries a standard
error and an interval, both read on the free scale its link defines and
mapped back
([`statmod_hyper_vcov()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)).
One chosen by
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) over
a kinked penalty is the argument of a minimum over a grid, so the row
names the criterion and leaves the remaining columns empty: there is no
curvature at such a point to read a standard error from. One the caller
set is marked fixed.

**Which test the statistic column reports is `test`'s answer.** By
default it is Wald's, the estimate over its standard error, which is the
only one of the four that needs no refit and is read off the
unrestricted fit alone. The other three – the likelihood ratio, Rao's
score and Terrell's gradient – hold the coefficient at the null value
and refit, one refit per row, and are reported as the SIGNED ROOT of
their \\\chi^2_1\\ statistic with the sign of \\\hat\beta_j - b\\, so
that the column carries one kind of quantity whatever was asked for and
Wald's own \\z\\ is the case where that root is exact. The p-value is
the \\\chi^2_1\\ one either way, and for Wald the two readings agree
identically.

Only the statistic and its p-value move with `test`: the interval stays
the Wald one, an inverted interval being some six times the cost of a
statistic per row.
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
takes `method` for that, one coefficient at a time. A row the restricted
fit is not defined for – a coefficient under a KINKED penalty, an
aliased one, any coefficient of a model carrying a structural term –
reports `NA` rather than falling back on Wald's, so the column never
carries two tests at once. Every other row is tested, a smooth's own
coordinates and a random effect's among them, on the penalized objective
and with the reading
[`statmod_stat_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_stat_at.md)
states.

**What a Wald p-value means here depends on the row**, and the summary
says which is which, in place of printing one column and leaving it at
that. For an unpenalized coefficient it is the usual thing. For a
coefficient in a penalized block it is conditional on the smoothing
parameter, which was not estimated jointly with it, and it does not
account for the shrinkage of the estimate towards zero. For a
coefficient a kinked penalty selected, the row exists only because that
coefficient survived the selection, and a naive interval there
under-covers.

**The degrees of freedom** are the effective ones, summed over the
terms, so a penalized term counts what it spends instead of how many
columns it has. The information criteria are built on that count.

## Reading the foot

Below the criteria the summary describes the point the fit stopped at
with up to five numbers from
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md).
No verdict is printed: the numbers are the point itself, and each is
dimensionless, so it does not move with the units of a coefficient or
with the sample size.

- `inner max |grad|/se`:

  the largest \\\|g_j\|/\sqrt{A\_{jj}}\\ over the coefficients, where
  \\g\\ is the gradient of the penalized objective and \\A\\ its
  curvature, the penalized information. \\1/\sqrt{A\_{jj}}\\ is the
  standard error of coefficient \\j\\ CONDITIONAL on the others, which
  is smaller than the marginal one the tables print, and the ratio
  equals \\\sqrt{2\Delta_j}\\ with \\\Delta_j = g_j^2/(2A\_{jj})\\ the
  log-likelihood a Newton step in that coordinate alone would still buy.
  At a mode it is small: fits located to their optimum read between
  1e-11 and 1e-2, a fit stopped after one iteration 0.84. For a model
  carrying a score-driven filter it covers the filter's own parameters
  as well. It leaves out a coefficient a kinked penalty set to zero, a
  frozen working block of
  [`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
  or
  [`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
  and an aliased column.

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
  leaving zero. A reading of 0.89 says the coefficient the data pull
  hardest is held by the kink with 11 per cent to spare, so with the
  other coefficients fixed it would leave zero at a \\\lambda\\ about 11
  per cent smaller; along a path the others move too, and the threshold
  moves with them. It is not a test of significance, only of whether the
  zeros agree with the \\\lambda\\ of the fit, and values close to one
  are ordinary where \\\lambda\\ was chosen on a grid, the next
  coefficient to enter being close to its threshold by construction.

- `outer max |grad|/se` and `min eigen`:

  the same two readings for the criterion that estimated the
  hyperparameters, on the free scale their links define, printed only
  where
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  or
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  ran. The curvature is the negated Hessian of the criterion. A
  hyperparameter that ran to the edge of its chart is left out, its
  curvature collapsing with its gradient; where every one did the line
  reads `not available: at a boundary`. A hyperparameter chosen along a
  path has no gradient, and there is no outer line.

## See also

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md),
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md),
which asks the same four of ONE coefficient against a value that need
not be zero, and
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
#>   (Intercept)     1.063 0.07514 14.14 < 1e-16 0.9156 1.210
#>   x               1.903 0.13020 14.62 < 1e-16 1.6480 2.158
#> 
#> === sigma   [log link]
#> 
#> Parametric terms
#>                estimate     se        z         p   lower   upper
#>   (Intercept)   -0.9647 0.1386 -6.96100 3.375e-12 -1.2360 -0.6931
#>   x              0.0216 0.2389  0.09041     0.928 -0.4466  0.4898
#> 
#> 95% intervals, bayesian variance
#> conditional log-likelihood -55.844138    effective df 4.00
#> cAIC 119.688    cBIC 130.838
#> fitted in 34 ms
#> inner   max |grad|/se 3.1e-06   min eigen 0.11
```
