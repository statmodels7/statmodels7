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
  structural = NULL,
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

  Optional named list, one entry per distribution parameter, each a
  numeric vector, a single number or a function of the coefficient
  count. See the details. `NULL` draws every coefficient.

- structural:

  Optional named list of a structural term's own parameters, on the
  scale
  [`modelterms7::term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.html)
  names. Strongly recommended when the formula carries such a term.

- sd:

  The standard deviation of the drawn coefficients, `1` by default. Read
  only for the coefficients `par` does not fix.

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

## The coefficients

`par = NULL` draws every one from `rnorm(1, 0, sd)`, which on the link
scale gives predictors of order one. A named list fixes them instead,
one entry per distribution parameter, and an entry may be:

- a numeric vector, as long as that equation has coefficients;

- a single number, used for every coefficient of the equation;

- a **function** of the coefficient count, called once and returning
  that many values.

The function is how a structured truth is written without a vocabulary
for it. `function(k) rnorm(k, 0, 0.3)` is a random effect with its own
standard deviation, and `function(k) c(1.5, -2, rep(0, k - 2))` is a
sparse truth for a lasso to find. A function answering with the wrong
count is refused, R being willing to recycle it into a different model.
A parameter left out of the list is drawn.

## A term with state

Simulated through
[`modelterms7::term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.html),
so the recursion that generates is the term's own. A score-driven term
draws the response as it runs, its level at one time driven by the score
at the time before; a latent chain draws its path from the stationary
law the likelihood is written with; a marginal break-point term draws
each group's positions from their prior. What each drew comes back in
`latent`, and that is what a recovery check compares against.

Such a term's own parameters are not coefficients of any equation, so
they are named through `structural`, never through `par`, on the scale
[`modelterms7::term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.html)
names: a loading is the loading, not its logarithm, a persistence is the
partial autocorrelation the chart carries. A formula holds at most one
such term, so no key is needed.

Left unnamed they take the term's own starting values, which are
deliberately weak. A score-driven term starts at a loading near 0.1, and
the series then has almost no dynamics: measured, its level ranged over
0.64 against 2.40 at named parameters. Name them, or the simulation is
of a model close to the one with no term at all.

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
#> [1] 1.255760 3.065784 6.383685

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
#> [1] 2.026176 2.036889 2.078646 2.068472 2.041046

# a score-driven series, its own parameters named
sim5 <- rstatmod(y ~ 0 + gas(p = 1, q = 1, time = t),
                 distributions7::gaussian1_distrib(),
                 data.frame(t = 1:100), par = list(sigma = 0),
                 structural = list(omega = 0.4, alpha1 = 0.3,
                                   pacf1 = 0.6))
head(sim5$latent, 3)
#> [1] 1.0000000 0.9164683 0.7841508
```
