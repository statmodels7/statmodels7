# The Coefficients of a Fitted Model

One named vector per distribution parameter: the quantities the model is
written in, or the coordinates it was estimated on.

## Usage

``` r
# S3 method for class 'StatmodFit'
coef(object, readable = TRUE, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- readable:

  `TRUE`, the default, reports the quantities the model is written in.
  `FALSE` reports the coordinates it was estimated on.

- ...:

  Unused.

## Value

A named list with one entry per distribution parameter, in the family's
order, each a named numeric vector. The names are the coefficient
labels, composed from each term's own label, so two terms of one kind in
one formula stay apart.

## What `readable` moves, and what it does not

Most coefficients are the same under either reading. A coefficient of a
linear predictor is what it is. Two kinds of parameter are reported
under a different name from the one they are carried under, and those
are the ones this argument moves.

A **break-point** term is fitted through a working pair and its position
is read off that pair. A discontinuous term carries \\g\\ and reports
\\\psi = -g/\delta\\, so at `readable = FALSE` the vector holds a number
that is no quantity of the model at all.

A **score-driven** term's persistence rides a partial autocorrelation,
the stationary region not being a box, and what the literature calls
\\\beta_j\\ is the autoregressive coefficient the whole chart produces.
At \\q = 2\\ a fit reporting \\\beta_1 = 0.761\\ has a free coordinate
of \\\mathrm{pacf}\_1 = 0.857\\.

Where a term declares no quantities of its own the coordinates stand. So
do the coefficients of a parameter developed over covariates: a
development is a vector with no single value to report, so the term
declares nothing for it.

## A structural term is present under both readings

It contributes no design columns, so its parameters are in no block.
They used to be in neither reading, and a model whose whole predictor is
a score-driven filter answered `numeric(0)`.

## When to ask for `FALSE`

`readable = FALSE` gives the vector the fit was estimated on, in the
order and under the names
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
is indexed by, with a structural term's part on the unconstrained scale
its charts define. That is what a caller feeding a fit back to
[`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
or to an optimizer needs.

Hyperparameters are not coefficients and are not here.
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
reports them.

## See also

[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
for the hyperparameters,
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the variance matrix indexed by the `readable = FALSE` names,
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
for the two together,
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html)
for what a term declares.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(200, 0, 10)))
d$y <- 0.3 * d$x + 1.5 * pmax(d$x - 6, 0) + rnorm(200, 0, 0.4)
fit <- statmod(y ~ modelterms7::seg(x, psi = 4),
               distributions7::gaussian1_distrib(), d)

# The break-point is a quantity of the model, near the true 6.
coef(fit)$mu
#> (Intercept)    seg.beta  seg.gamma1    seg.psi1 
#> -0.05693875  0.31903803  1.54956901  6.13040161 

# The coordinates the fit ran on are the same numbers here, seg() being
# continuous, and vcov() is indexed by these names.
coef(fit, readable = FALSE)$mu
#> (Intercept)    seg.beta  seg.gamma1    seg.psi1 
#> -0.05693875  0.31903803  1.54956901  6.13040161 
rownames(vcov(fit))
#> [1] "mu:(Intercept)"    "mu:seg.beta"       "mu:seg.gamma1"    
#> [4] "mu:seg.psi1"       "sigma:(Intercept)"
```
