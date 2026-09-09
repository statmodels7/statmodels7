# Simulate a Response From a Written Model

Takes a formula and a distribution, draws coefficients or uses the ones
given, and returns the data drawn from that model together with the
truth behind it.

## Usage

``` r
rstatmod(
  formula,
  distrib,
  data = NULL,
  n = NULL,
  n_sim = 1,
  par = NULL,
  sd = 1,
  offsets = NULL,
  covariates = NULL
)
```

## Arguments

- formula:

  The model formula, as
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  takes it, with the equations separated by `|`. Its left-hand side must
  be a symbol.

- distrib:

  A distributions7 distribution object, which decides how many equations
  there are and what is drawn from.

- data:

  A data frame of covariates, or `NULL` where `n` is given.

- n:

  The number of observations, where `data` is `NULL`. One of the two is
  required.

- n_sim:

  How many data sets to draw. `1` by default.

- par:

  Optional named list of what to hold, over the whole model: a
  distribution parameter, one of its coefficients or a group of them, a
  structural term's own parameter or a group of those. Each entry is a
  numeric vector, a single number or a function of the count. See the
  details. `NULL`, the default, draws everything.

- sd:

  The width of the draws, `1` by default. A quantity that rides a chart
  – a structural parameter, a prior's own scale – is drawn at half of
  it, the chart carrying a width of one onto most of its parameter's
  range.

- offsets:

  Optional named list of offsets, one per parameter, as
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  takes them.

- covariates:

  Optional named list of functions of the observation count, one per
  covariate, drawn afresh at every replicate. A drawn factor must have
  its levels fixed.

## Value

An object of class `"StatmodSim"`, a list of seven:

- `data`:

  the data frame with the response column added, named after the
  formula's left-hand side.

- `par`:

  the coefficients used, drawn or given, named as the design names them.

- `theta`:

  the distribution's parameters at every observation, a named list with
  one vector per parameter.

- `latent`:

  what a term with state drew, or `NULL`.

- `structural`:

  such a term's own parameters, or `NULL`.

- `hyper`:

  a data frame of the hyperparameters, one row per penalty and name,
  with the value drawn or held. Its columns are those of
  [`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
  as far as they mean the same thing, so a study compares the two
  directly.

- `n_sim`:

  as supplied.

- `call`:

  the matched call.

With `n_sim > 1` the fields `data`, `theta` and `latent` are lists of
that length.

## What it is for

Data whose truth is known: write the model, draw from it, fit it back,
and see whether the fit recovers what was put in. That is the shape of a
simulation study, and of a check on a term one has just written.

This is not
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html), which
draws from a model already fitted. The `r` prefix is R's own for a
random draw, so the two names stay apart.

A covariate needs no declaring. A factor becomes its contrasts and a
numeric stays itself, the design coming from the same interpreter a fit
uses.

## The truth comes back beside the data

A simulation study compares against the coefficients, the parameters
they gave and whatever a term with state drew, so the result is a list
holding all of them. Pass `sim$data` where a data frame is wanted.

They were attributes of the data frame until version 0.88.0, and that
was worse than it looks: an attribute survives a row subset without
being subset itself, so `sim[1:10, ]` silently kept a `theta` of the
original length, while [`subset()`](https://rdrr.io/r/base/subset.html)
and [`merge()`](https://rdrr.io/r/base/merge.html) dropped it
altogether.

## The predictor is assembled as a fit assembles it

Through
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md),
so the simulated data come from the model that was written, never from a
linearization of it.

A term whose block moves with its coefficients,
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
or
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
contributes `term_value()` at the coefficients supplied, not its block
times them. The two differ by the whole nonlinearity, and an earlier
version of this function used the second: measured, the gap is 3 on a
`seg()`, 4.05 on an `nl()` and a missing value on a `jseg()`.

## Data, or a row count

`data` carries the covariates. A model with none, a pure time series
say, needs only `n`. One of the two must be given, and where both are
given they must agree.

## Fixed covariates or drawn ones

`covariates` takes one function of the observation count per column, as
`par` takes one per equation, and they are drawn afresh at every
replicate.

The choice decides what a study measures. With `data` the covariates are
the same throughout, so what is measured is the estimator's behavior
**conditional** on that design; with `covariates` it is measured over
the design as well. Neither is more correct, and a study should say
which it ran. Measured on a simple regression at \\n = 40\\, the slope's
standard deviation is 0.1465 under the first and 0.1443 under the
second, so the distinction is about what a study claims, never about a
large numerical difference.

A drawn factor is refused unless its levels are fixed. Drawn freely it
loses a level on some replicate, and the coefficients drawn against the
first design would then be recycled into a different model.

## Several replicates

`n_sim` draws that many data sets. The truth is drawn **once** and
shared: a study over replicates measures the variability of an estimator
at a set of parameters, so the replicates differ in what is random,
never in what is being estimated. Varying the truth as well is a loop
over calls and reads differently.

With `n_sim > 1` the per-replicate fields, `data`, `theta` and `latent`,
come back as lists of that length, while `par` and `structural` stay
single.

## Everything is drawn

A call that names nothing draws the whole truth, so writing the model is
the whole of what getting data from it takes. The rule is that whoever
knows what a quantity means draws it:

- a coefficient of a design column has no other owner and comes from
  `rnorm(1, 0, sd)`, which on the link scale gives predictors of order
  one;

- a coordinate some penalty covers is drawn from that penalty read as a
  prior, through
  [`penalties7::penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.html).
  A Gaussian random effect gives Gaussian effects, a lasso Laplace ones
  and a heavy-tailed prior heavy-tailed ones. The prior's own scale is
  drawn as well, on the chart the penalty carries for it, and comes back
  in `hyper`;

- a structural term's own parameters are drawn by the term, through
  [`modelterms7::term_draw()`](https://statmodels7.github.io/modelterms7/reference/term_draw.html),
  which knows the chart each one rides. A loading stays positive and a
  persistence stationary whatever comes out;

- a coefficient of a design column whose meaning only the term knows is
  drawn by the term, through
  [`modelterms7::term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.html),
  last of all so that what it writes is what survives. A break-point is
  the case: it is a position on the covariate's own axis, so a normal of
  width `sd` lands outside the data as often as not and the confinement
  then pins it to the interval's edge, where one of the two segments
  holds a twentieth of the rows. Measured over fifty groups with the
  covariate uniform on \\(0, 1)\\, thirty-eight of the fifty came back
  pinned and twelve strictly interior.

A hyperparameter the term holds is used rather than drawn, so
`s(x, bspline_smooth(), hyper = c(lambda = 2))` simulates at the
smoothing it names. A prior whose coordinates a term drew is reported at
the width the term used, the values there no longer being that prior's
own draw – and drawn Gaussian at that width, so a heavy-tailed prior
over a break-point keeps its family for the fitting and not for the
simulation.

What no prior reaches falls back to the plain draw, and it is a short
list: SCAD and MCP are improper by construction, an anisotropic tensor
smooth is flat along its null space, and a covariance class spanning a
filter and an equation at once is in neither vector on its own.

## Holding what you care about

`par` is one named list over the whole model, and a key may be

- a distribution parameter, `mu`, which is that whole equation;

- one of its coefficients or a group of them, `mu.(Intercept)` or
  `mu.random`, a group being a name the members extend at a dot;

- a structural term's own parameter, `alpha1`, or a group of those,
  `omega.random`.

A value is a vector of that key's own length, a single number used for
all of them, or a **function** of the count. The function is how a
structured truth is written without a vocabulary of its own:
`function(k) rnorm(k, 0, 0.4)` is a random effect at a scale one chose
and `function(k) c(1.5, -2, rep(0, k - 2))` is a sparse truth for a
lasso to find. A function answering with the wrong count is refused, R
being willing to recycle it into a different model, and a key that
reaches nothing is refused with what the model does carry.

A structural parameter is named on the scale a reader knows, which is
[`modelterms7::term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.html)'s:
a loading is the loading and not its logarithm, a persistence the
partial autocorrelation its chart carries. A parameter a subformula
DEVELOPS is different, and it has to be: its coordinates are the
coefficients of that development, which act on the unconstrained scale
of the parameter's own chart, so `alpha1` is a loading and
`alpha1.random.3` is a group's departure on the log scale that loading
rides. That is what keeps every group's loading positive whatever the
departure is.

## A term with state

Simulated through
[`modelterms7::term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.html),
so the recursion that generates is the term's own. A score-driven term
draws the response as it runs, its level at one time driven by the score
at the time before; a latent chain draws its path from the stationary
law the likelihood is written with; a marginal break-point term draws
each group's positions from their prior. What each drew comes back in
`latent`, and that is what a recovery check compares against.

Such a term's own parameters are not coefficients of any equation and
are named in `par` alongside them, a formula holding at most one such
term so that no key is needed. They are drawn like everything else,
which they were not until version 0.111.0: a filter whose level was
developed over a hundred groups took its starting values, and those are
zero for every deviation, so the panel came out with no heterogeneity
between groups at all – three distinct values of the mean over six
hundred observations.

## The response's name

The formula's left-hand side, which must be a symbol. `log(y) ~ x` is
refused: the model generates values of `log(y)`, and no column could
honestly be called either name. A censored response is refused too, for
the reason
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
refuses one.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
to fit what this draws,
[`simulate.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/simulate.StatmodFit.md)
to draw from a model already fitted,
[`print.StatmodSim()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodSim.md)
for the printed form.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(50), g = factor(rep(c("a", "b"), 25)))

# coefficients drawn
sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
sim$par
#> $mu
#>  (Intercept)            x           gb 
#>  0.291446236 -0.443291873  0.001105352 
#> 
#> $sigma
#> (Intercept) 
#>  0.07434132 
#> 

# or given, and recovered by a fit
sim2 <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                 par = list(mu = c(1, 2), sigma = log(0.3)))
coef(statmod(y ~ x, distributions7::gaussian1_distrib(), sim2$data))
#> $mu
#> (Intercept)           x 
#>    1.030843    1.935854 
#> 
#> $sigma
#> (Intercept) 
#>   -1.165592 
#> 

# a sparse truth, written as a function of the coefficient count
dd2 <- as.data.frame(matrix(rnorm(50 * 6), 50, 6))
sim3 <- rstatmod(y ~ lasso(~ V1 + V2 + V3 + V4 + V5 + V6),
                 distributions7::gaussian1_distrib(), dd2,
                 par = list(mu = function(k) c(2, -1.5, rep(0, k - 2)),
                            sigma = log(0.3)))
head(sim3$data$y, 3)
#> [1] 2.108603 3.625229 6.596983

# a model with no covariates at all
sim4 <- rstatmod(y ~ 1, distributions7::gaussian1_distrib(), n = 20)
nrow(sim4$data)
#> [1] 20

# a study over replicates, the covariates drawn afresh at each one
study <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), n = 80,
                  n_sim = 5, par = list(mu = c(1, 2), sigma = log(0.5)),
                  covariates = list(x = function(n) runif(n, -2, 2)))
length(study$data)
#> [1] 5
vapply(study$data, function(d) coef(statmod(
  y ~ x, distributions7::gaussian1_distrib(), d),
  readable = FALSE)$mu[[2L]], numeric(1))
#> [1] 2.017968 2.036897 2.013920 2.043438 1.935684

# a score-driven series, its own parameters named beside the equations'
sim5 <- rstatmod(y ~ 0 + gas(p = 1, q = 1, time = t),
                 distributions7::gaussian1_distrib(),
                 data.frame(t = 1:100),
                 par = list(sigma = 0, omega = 0.4, alpha1 = 0.3,
                            pacf1 = 0.6))
head(sim5$latent, 3)
#> [1] 1.0000000 0.6960904 0.6137012

# and a panel whose level varies by group: naming nothing draws the
# deviations from the prior random() declares, at a scale drawn too
pan <- data.frame(id = factor(rep(1:6, each = 10)), t = rep(1:10, 6))
sim6 <- rstatmod(y ~ 0 + gas(p = 1, q = 1, omega ~ 1 + random(~1 | id),
                             by = id, time = t),
                 distributions7::gaussian1_distrib(), pan)
sim6$hyper
#>   parameter                                                              term
#> 1        mu gas(p = 1, q = 1, omega ~ 1 + random(~1 | id), by = id, time = t)
#>    name     value  held
#> 1 sigma 0.6269308 FALSE
length(unique(round(sim6$theta$mu, 6))) > 6
#> [1] TRUE
```
