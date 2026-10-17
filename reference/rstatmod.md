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

  A distributions7 distribution object, which determines how many
  equations there are and the distribution from which the data are
  drawn.

- data:

  A data frame of covariates, or `NULL` where `n` is given.

- n:

  The number of observations, a positive whole number. It is required
  where `data` is `NULL`, and where `data` is given it must equal the
  number of rows.

- n_sim:

  How many data sets to draw, a positive whole number. The default is
  `1`.

- par:

  Optional named list of what to hold, over the whole model: a
  distribution parameter, one of its coefficients or a group of them, a
  structural term's own parameter or a group of those. Each entry is a
  numeric vector, a single number or a function of the count. See the
  details. `NULL`, the default, draws everything.

- sd:

  The width of the draws, `1` by default. A quantity that rides a chart
  (a structural parameter, a prior's own scale) is drawn at half of it,
  because a chart carries a width of one onto most of its parameter's
  range.

- offsets:

  Optional named list of offsets, one per parameter, as
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  takes them.

- covariates:

  Optional named list of functions of the observation count, one per
  covariate, drawn afresh at every replicate. A drawn factor should have
  its levels fixed, since a replicate that loses a level stops the
  simulation with an error.

## Value

An object of class `"StatmodSim"`, a list with the elements:

- `data`:

  the data frame, to which the response column has been added under the
  name of the formula's left-hand side.

- `par`:

  the coefficients used, drawn or given, named as the design names them.

- `theta`:

  the distribution's parameters at every observation, a named list with
  one vector per parameter.

- `latent`:

  what a term with state drew (the filtered level of a score-driven
  term, the states of a latent chain, or the positions of a marginal
  break-point term), or `NULL`.

- `structural`:

  such a term's own parameters, or `NULL`.

- `hyper`:

  a data frame of the hyperparameters, with one row per penalty and
  name; each row gives the value drawn or held. Its columns are those of
  [`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
  as far as they mean the same thing, and `value` takes the place of
  `estimate`, so a study compares the two directly.

- `n_sim`:

  the number of data sets, as an integer.

- `call`:

  the matched call.

With `n_sim > 1` the fields `data` and `theta` are lists of that length,
and `latent` is such a list when the model has a term with state and
`NULL` otherwise; `par`, `structural` and `hyper` stay single.

## What it is for

Simulated data have a known truth: write the model, draw from it, fit it
back, and see whether the fit recovers what was put in. That is the
shape of a simulation study, and of a check on a term that has just been
written.

This is not
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html), which
draws from a model already fitted. The `r` prefix is the one that R uses
for a random draw, so the two names stay apart.

A covariate needs no declaration. A factor becomes its contrasts and a
numeric variable stays as it is, and the design is built by the
interpreter that a fit uses.

## The returned truth

A simulation study compares against the coefficients, the parameters
that they gave and whatever a term with state drew, so the result is a
list holding all of them. Pass `sim$data` where a data frame is wanted.

The truth is kept as elements of a list and not as attributes of the
data frame, because an attribute is not subset with the rows:
`sim[1:10, ]` would keep a `theta` of the original length, while
[`subset()`](https://rdrr.io/r/base/subset.html) and
[`merge()`](https://rdrr.io/r/base/merge.html) drop it.

## Assembly of the predictor

The predictor is built through
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md),
so the simulated data come from the model that was written and not from
a linearization of it.

A term whose block moves with its coefficients,
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
or
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
contributes `term_value()` at the coefficients supplied, not its block
times them. The two differ by the whole nonlinearity, and the simulation
uses the first.

## Data, or a row count

`data` carries the covariates. A model with none, a pure time series
say, needs only `n`. One of the two must be given, and where both are
given they must agree.

## Fixed covariates or drawn ones

`covariates` takes one function of the observation count per column, as
an entry of `par` may be, and the columns are drawn afresh at every
replicate. A column of `data` under the same name is overwritten, and
the columns of `data` that are not drawn are the same in every
replicate.

The choice determines what a study measures. With `data` the covariates
are the same throughout, so the study measures the behavior of the
estimator conditional on that design; with `covariates` it measures that
behavior over the design as well. The two choices are equally correct,
and a study should state which one it ran.

A drawn factor should have its levels fixed. A factor drawn with free
levels can lose a level at some replicate. The design then has another
shape than the one for which the coefficients were drawn, and the
simulation stops at that replicate with an error.

## Several replicates

`n_sim` draws that many data sets. The truth is drawn once and shared by
the replicates: a study over replicates measures the variability of an
estimator at a given set of parameters, so the replicates differ in what
is random and not in what is being estimated. Varying the truth as well
requires a loop over calls and reads differently.

With `n_sim > 1` the per-replicate fields `data` and `theta` come back
as lists of that length, and `latent` does so when the model has a term
with state (otherwise it is `NULL`), while `par`, `structural` and
`hyper` stay single.

## The quantities that are drawn

If the user names no quantity, the call draws the whole truth, so
writing the model is enough to get data from it. Each quantity is drawn
by the component that defines its meaning:

- a coefficient of a design column has no other owner and comes from
  `rnorm(1, 0, sd)`, which on the link scale gives predictors of order
  one;

- a coordinate that some penalty covers is drawn from that penalty read
  as a prior, through
  [`penalties7::penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.html).
  A Gaussian random effect gives Gaussian effects, a lasso Laplace ones
  and a heavy-tailed prior heavy-tailed ones. The prior's own scale is
  drawn as well, on the chart that the penalty carries for it, and comes
  back in `hyper`;

- a structural term's own parameters are drawn by the term, through
  [`modelterms7::term_draw()`](https://statmodels7.github.io/modelterms7/reference/term_draw.html),
  on the chart that each one rides. A loading stays positive and a
  persistence stationary whatever comes out;

- a coefficient of a design column whose meaning only the term defines
  is drawn by the term, through
  [`modelterms7::term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.html),
  last of all, so that the values that it writes survive. A break-point
  is such a coordinate: it is a position on the covariate's own axis, so
  a normal of width `sd` often lands outside the data, and the term then
  confines it, which pins it to an edge of the interval, where one of
  the two segments holds a small share of the rows.

A hyperparameter that the term holds is used and not drawn, so
`s(x, bspline_smooth(), hyper = c(lambda = 2))` simulates at the
smoothing parameter that it names. If a term draws the coordinates that
a prior covers, the prior is reported at the width that the term used,
since those values are no longer the prior's own draw. They are drawn
Gaussian at that width, so a heavy-tailed prior over a break-point keeps
its family for the fitting but not for the simulation.

A coordinate that a prior does not reach falls back to the plain draw.
The cases are SCAD and MCP, which are improper by construction; an
anisotropic tensor smooth, which is flat along its null space; a
covariance class spanning a filter and an equation at once, which is in
neither vector on its own; and a coordinate for which the draw of the
penalty is missing.

## Fixing quantities by name

`par` is one named list over the whole model, and a key may be

- a distribution parameter, `mu`, which is that whole equation;

- one of its coefficients or a group of them, `mu.(Intercept)` or
  `mu.random`, where a group is a name that the members extend at a dot;

- a structural term's own parameter, `kappa1`, or a group of those,
  `omega.random`.

A value is a vector of that key's own length, a single number used for
all of them, or a function of the count. A function writes a structured
truth without a vocabulary of its own: `function(k) rnorm(k, 0, 0.4)` is
a random effect at a scale chosen by the caller and
`function(k) c(1.5, -2, rep(0, k - 2))` is a sparse truth for a lasso to
find. A function that returns the wrong number of values or a missing
value is rejected, because R would recycle a vector of the wrong length
into a different model. A key that addresses nothing is rejected with a
message that lists what the model does carry.

A structural parameter is named on the scale of
[`modelterms7::term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.html):
a loading is the loading and not its logarithm, and a persistence is the
partial autocorrelation that its chart carries. A parameter that a
subformula develops is different. Its coordinates are the coefficients
of that development, which act on the unconstrained scale of the
parameter's own chart, so `kappa1` is a loading and `kappa1.random.3` is
a group's departure on the log scale that the loading rides. This keeps
every group's loading positive whatever the departure is.

## A term with state

The simulation runs through
[`modelterms7::term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.html),
so the recursion that generates the data is the term's own. A
score-driven term draws the response as it runs, and its level at one
time is driven by the score at the time before; a latent chain draws its
path from the stationary law with which the likelihood is written; a
marginal break-point term draws each group's positions from their prior.
The values that each term drew come back in `latent`, and a recovery
check compares a fit with them.

The parameters of such a term are not coefficients of any equation. They
are named in `par` beside the equations' coefficients, and since a
formula holds at most one such term, the names carry no term key. They
are drawn like every other coefficient.

## The response's name

The left-hand side of the formula must be a symbol. `log(y) ~ x` is
rejected because the model generates values of `log(y)`, which cannot be
stored in a column under either name. A censored response is rejected
too, for the reason that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
rejects one: there is no censored likelihood.

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

# The coefficients are drawn.
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

# The coefficients are given, and a fit recovers them.
sim2 <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                 par = list(mu = c(1, 2), sigma = log(0.3)))
coef(statmod(y ~ x, distributions7::gaussian1_distrib(), sim2$data))
#> $mu
#> (Intercept)           x 
#>    1.030843    1.935854 
#> 
#> $sigma
#> (Intercept) 
#>   -1.145181 
#> 

# A sparse truth is written as a function of the coefficient count.
dd2 <- as.data.frame(matrix(rnorm(50 * 6), 50, 6))
sim3 <- rstatmod(y ~ lasso(~ V1 + V2 + V3 + V4 + V5 + V6),
                 distributions7::gaussian1_distrib(), dd2,
                 par = list(mu = function(k) c(2, -1.5, rep(0, k - 2)),
                            sigma = log(0.3)))
head(sim3$data$y, 3)
#> [1] 2.108603 3.625229 6.596983

# This model has no covariates at all.
sim4 <- rstatmod(y ~ 1, distributions7::gaussian1_distrib(), n = 20)
nrow(sim4$data)
#> [1] 20

# A study over replicates draws the covariates afresh at each one.
study <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), n = 80,
                  n_sim = 5, par = list(mu = c(1, 2), sigma = log(0.5)),
                  covariates = list(x = function(n) runif(n, -2, 2)))
length(study$data)
#> [1] 5
vapply(study$data, function(d) coef(statmod(
  y ~ x, distributions7::gaussian1_distrib(), d),
  readable = FALSE)$mu[[2L]], numeric(1))
#> [1] 2.017968 2.036897 2.013920 2.043438 1.935684

# In a score-driven series, the term's own parameters are named beside the
# equations'.
sim5 <- rstatmod(y ~ 0 + gas(p = 1, q = 1, time = t),
                 distributions7::gaussian1_distrib(),
                 data.frame(t = 1:100),
                 par = list(sigma = 0, omega = 0.4, kappa1 = 0.3,
                            pacf1 = 0.6))
head(sim5$latent, 3)
#> [1] 1.0000000 0.6960904 0.6137012

# In a panel whose level varies by group, naming nothing draws the
# deviations from the prior that random() declares, at a scale drawn too.
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
