# The Profile Likelihood of Coefficients

Returns, for each combination of values of one or several coefficients,
the maximum of the penalized log-likelihood over all the other
coefficients when those are held at the values. The result is a data
frame, one row per combination, which serves a plot, an interval or a
diagnostic; the function computes none of those.

## Usage

``` r
statmod_profile(
  fit,
  which,
  values = profile_range(width = 3, n = 21),
  method = profile_grid(),
  hyper = c("held", "refit"),
  cores = 1
)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- which:

  A character vector of unique labels `"parameter:coefficient"`, of one
  equation or of several.

- values:

  Either one
  [`profile_range()`](https://statmodels7.github.io/statmodels7/reference/profile_range.md),
  which applies to every coefficient of `which`, or a list. The named
  elements of a list carry names of `which` and are numeric vectors (the
  values, on the scale described above) or
  [`profile_range()`](https://statmodels7.github.io/statmodels7/reference/profile_range.md)
  objects. At most one element without a name sets the default for the
  coefficients that no element names, and without it that default is
  [`profile_range()`](https://statmodels7.github.io/statmodels7/reference/profile_range.md).
  Two unnamed elements, or a name that is not in `which`, are errors.
  The values of each coefficient are sorted increasingly and repeated
  values are dropped. With several coefficients the table has one row
  for each combination of values, and the first coefficient varies
  fastest.

- method:

  A
  [profile_method](https://statmodels7.github.io/statmodels7/reference/profile_method.md):
  [`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md).

- hyper:

  `"held"` or `"refit"`: whether the hyperparameters stay at the values
  of the fit or are searched again at every point.

- cores:

  The number of separate R processes over which the refits are
  distributed, a single whole number (default `1`, no process besides
  this one). It is distinct from the `threads` of
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
  which counts threads inside one process; a refit that runs in a
  separate process is sequential. Under `pkgload` the run stays
  sequential.

## Value

A data frame (`check.names = FALSE`) with one column for each
coefficient of `which`, named by its label, then `profile_loglik`,
`profile_penalty`, `profile_pen_loglik`, `outer_criterion` and
`converged`. Its attributes are `maximum`, a one-row data frame with the
same columns at the unrestricted fit (the numbers of the fit itself,
which under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
can lie slightly below the nearest rows of the table, because at a row
the criterion does not integrate the held coefficients), `method`, the
[profile_method](https://statmodels7.github.io/statmodels7/reference/profile_method.md)
object, `hyper`, and `criterion`, the outer criterion object of the fit
(`NULL` if it has none). The hyperparameters and the other coefficients
of each point are not returned.

## What the columns report

The columns `profile_loglik`, `profile_penalty` and `profile_pen_loglik`
belong to the INNER criterion, the penalized log-likelihood of the
coefficients at the hyperparameters in force. At the point where the
coefficients are held at the values of the row, `profile_loglik` is the
log-likelihood without the penalty, `profile_penalty` is the penalty
with the sign that makes the two add up, that is the logarithm of the
prior density of the coefficients, \\-\tfrac12 P\_\lambda(\beta)\\ with
its normalizing constant, and `profile_pen_loglik` is their sum, which
is the quantity that the refit maximizes. The constant depends on the
hyperparameters and not on the coefficients, so the penalty is not
necessarily negative, and with the hyperparameters held its changes
along the grid are those of the quadratic form. A model without penalty
has `profile_penalty` equal to zero, and the two log-likelihoods then
coincide.

The column `outer_criterion` belongs to the OUTER criterion
([`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)),
a different function of a different object: the value, on the sign and
the scale that the fit reports in `fit@criterion`, of the criterion that
chooses the hyperparameters or the coefficients that the criterion
estimates. It is not the profile of that criterion over the held
coefficients in any sense that the first three columns have for the
inner one, and it is `NA` when the point involved no outer search: when
the fit has no outer criterion, or when the hyperparameters are held
(`hyper = "held"`) and the criterion estimates no coefficient.

The column `converged` states that the point is at the maximum of its
refit. It is `TRUE` when the stopping rule of the refit fired, or when
the Newton decrement at the point is below the tolerance of
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
and in both cases the final inner fit converged. The second condition is
needed because the line search of the outer optimizer can exhaust its
backtracks at a point that is already the optimum, without recording
that its stopping rule fired. The row of the attribute `maximum` reports
instead `fit@converged`, which is only the flag of the stopping rule. A
refit that signals an error gives a row of `NA` with
`converged = FALSE`, and the function warns.

## The coefficients to hold: `which`

A label has the form `"parameter:coefficient"`: the distribution
parameter whose equation carries the coefficient, a colon, and the name
of the coefficient as the second column of the table of
[`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md)
gives it, which is the name that
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) print. The labels of
a fit are the row names of `summary(fit)`, for example
`"mu:(Intercept)"`, `"mu:x"` or `"sigma:s(z).1"`. Labels must be unique.

The coefficient of the effect of a covariate `z` on the parameter
`theta` of a
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
term, written `nl(~ a * exp(-theta * x), theta ~ 1 + z)` in the equation
of `mu`, has the label `"mu:nl.theta.z"`. The name after the colon is
`nl.`, the name of the parameter of the term, a dot and the name of the
coefficient of its subformula, so that the intercept of `theta` is
`"mu:nl.theta.(Intercept)"`. A parameter of the term without a
subformula is a coefficient by itself, as in `"mu:nl.a"`.

The position of a break-point of
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
is the coefficient `psi1` (then `psi2`, and so on) of the term, as in
`"mu:seg.psi1"`. When the position has a subformula, as in
`seg(w, psi ~ x)`, its coefficients are `"mu:seg.psi1.(Intercept)"` and
`"mu:seg.psi1.x"`. The names of all the coefficients of a fit are the
row names of
`statmodels7:::coef_labels(fit@spec, statmodels7:::statmod_design(fit@spec))`.

When two or more positions of the same `seg()` term without subformulas
are in `which`, the positions of a model are ordered and distinct: the
combinations of the grid in which a position is not larger than the
position before it are not evaluated and are not in the table.

## The scale

The values are on the scale of the coefficient in the stacked vector of
the fit, which is unconstrained: a coefficient of a linear predictor is
on the scale of the link, and so are the coefficients of the predictor
of a parameter of a nonlinear term, the positions of the break-points
and the parameters of a covariance structure. Values on another scale
must be mapped to it by the caller. A monotone map relabels the grid and
leaves the columns of the profile unchanged.

## The hyperparameters, and the coefficients of the criterion

With `hyper = "held"` the hyperparameters (the smoothing parameters and
the scales of the random effects) stay at the values of the fit, and
with `hyper = "refit"` the outer search is run again at every point. In
both cases, when the outer criterion estimates coefficients
([`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimate the unpenalized coefficients of every parameter except the
first), those coefficients are estimated by the criterion at every point
and are not carried to the joint mode of the penalized likelihood, which
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md)
does. A held coefficient that the criterion would have estimated is
held, and the criterion estimates the rest. The cost of a point is then
that of an outer search over those coordinates, and it exceeds the cost
of a restricted fit.

The outer optimizer is the default one of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
and the stopping rules, the inner method and the information are those
of the fit.

## Coefficients that cannot be held yet

A coefficient under a penalty with a kink (a lasso, a SCAD or an MCP), a
coefficient of a model with a structural term
([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html),
the marginal break-point term) and any coefficient with
`hyper = "refit"` in a model whose kinked penalties are tuned by a path
signal an error. An aliased coefficient has no estimate to profile and
signals an error.

## See also

[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md),
[`profile_range()`](https://statmodels7.github.io/statmodels7/reference/profile_range.md),
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md)
for the likelihood ratio at one value.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)

# the profile of the slope, on five standard errors
prof <- statmod_profile(fit, which = "mu:x",
                        values = profile_range(width = 2, n = 9))
prof
#>       mu:x profile_loglik profile_penalty profile_pen_loglik outer_criterion
#> 1 1.615742      -23.56354               0          -23.56354       -25.70965
#> 2 1.700842      -22.70458               0          -22.70458       -24.86501
#> 3 1.785942      -22.07563               0          -22.07563       -24.24654
#> 4 1.871043      -21.69183               0          -21.69183       -23.86914
#> 5 1.956143      -21.56280               0          -21.56280       -23.74225
#> 6 2.041243      -21.69183               0          -21.69183       -23.86914
#> 7 2.126343      -22.07563               0          -22.07563       -24.24654
#> 8 2.211443      -22.70458               0          -22.70458       -24.86501
#> 9 2.296543      -23.56354               0          -23.56354       -25.70965
#>   converged
#> 1      TRUE
#> 2      TRUE
#> 3      TRUE
#> 4      TRUE
#> 5      TRUE
#> 6      TRUE
#> 7      TRUE
#> 8      TRUE
#> 9      TRUE

# the log-likelihood falls by z^2 / 2 at z standard errors
attr(prof, "maximum")
#>       mu:x profile_loglik profile_penalty profile_pen_loglik outer_criterion
#> 1 1.956143      -21.57563               0          -21.57563       -24.59838
#>   converged
#> 1      TRUE

# two coefficients, with values chosen by the caller for the second
statmod_profile(fit, which = c("mu:x", "mu:(Intercept)"),
                values = list("mu:x" = profile_range(width = 1, n = 3),
                              "mu:(Intercept)" = c(0.8, 1, 1.2)))
#>       mu:x mu:(Intercept) profile_loglik profile_penalty profile_pen_loglik
#> 1 1.785942            0.8      -43.60960               0          -43.60960
#> 2 1.956143            0.8      -35.87738               0          -35.87738
#> 3 2.126343            0.8      -29.40230               0          -29.40230
#> 4 1.785942            1.0      -27.66303               0          -27.66303
#> 5 1.956143            1.0      -22.79422               0          -22.79422
#> 6 2.126343            1.0      -22.13489               0          -22.13489
#> 7 1.785942            1.2      -22.49681               0          -22.49681
#> 8 1.956143            1.2      -25.44694               0          -25.44694
#> 9 2.126343            1.2      -31.78035               0          -31.78035
#>   outer_criterion converged
#> 1       -43.60960      TRUE
#> 2       -35.87738      TRUE
#> 3       -29.40230      TRUE
#> 4       -27.66303      TRUE
#> 5       -22.79422      TRUE
#> 6       -22.13489      TRUE
#> 7       -22.49681      TRUE
#> 8       -25.44694      TRUE
#> 9       -31.78035      TRUE
```
