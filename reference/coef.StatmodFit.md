# The Coefficients of a Fitted Model

One named vector per distribution parameter: the quantities in which the
model is written, or the coordinates on which it was estimated.

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

  `TRUE`, the default, reports the quantities in which the model is
  written. `FALSE` reports the coordinates on which it was estimated.

- ...:

  Unused.

## Value

A named list with one entry per distribution parameter, in the family's
order, each a named numeric vector. The names are the coefficient
labels, composed from each term's own label, so two terms of one kind in
one formula stay apart.

## What `readable` moves, and what it does not

Most coefficients are the same under either reading: the coefficient of
a linear predictor is reported as it is estimated. Three kinds of term
report quantities that differ from the coordinates they are carried on,
and those are the parameters that this argument moves.

A **break-point** term is fitted through a working block, and its
quantities are read from that block. A continuous term holds each
quantity as a coefficient, and the two readings coincide. A
discontinuous term carries \\g\\ and \\\delta\\; where its position is
not held, it reports \\\psi = -g/\delta\\, so that at `readable = FALSE`
the vector holds \\g\\, which is no quantity of the model. Where the
term holds its position, as a sharp `jump()` or `jseg()` does after the
fit has settled it, the reported position is the held one.

A **score-driven** term carries its level and its loadings \\\kappa_j\\
through their own links, the log for a loading, and its persistence on
partial autocorrelations, the stationary region not being a box. Each
partial autocorrelation is carried on the unconstrained scale of its
link, the inverse hyperbolic tangent by default. What the literature
calls \\\xi_j\\ is the autoregressive coefficient that the whole chart
produces. With one persistence, \\\xi_1\\ is the partial
autocorrelation. With two, \\\xi_1 = \mathrm{pacf}\_1 (1 -
\mathrm{pacf}\_2)\\ and \\\xi_2 = \mathrm{pacf}\_2\\.

A **regime** term is estimated on the additive log-ratios of each row of
its transition matrix, and it reports the transition probabilities
\\p\_{ij}\\ in their place. Its levels and gaps are reported through
their own links.

Where a term declares no quantities of its own, the coordinates stand.
So do the coefficients of a parameter developed over covariates: a
development is a vector with no single value to report, so the term
declares nothing for it.

## A structural term is present under both readings

A structural term contributes no design columns, so its parameters are
in no block of the design. Both readings append them to the vector of
the equation, so that a model whose whole predictor is a score-driven
filter returns its filter parameters.

## The coordinates of the fit

`readable = FALSE` gives the coordinates on which the fit was estimated,
with a structural term's part on the unconstrained scale that its charts
define. For a model without a structural term, this vector has the order
of the rows of `vcov(object, readable = FALSE)`, and it can be passed as
the `par` argument of
[`loglik()`](https://statmodels7.github.io/statmodels7/reference/loglik.md)
or to an optimizer. A structural level that is held appears in the
vector at its held value and has no row in that matrix.

Hyperparameters are not coefficients and are not here.
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
reports them.

## See also

[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
for the hyperparameters,
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the variance matrix, whose rows with `readable = FALSE` are the
unheld coordinates of `coef(object, readable = FALSE)`,
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
# continuous, and vcov() has one row for each of them.
coef(fit, readable = FALSE)$mu
#> (Intercept)    seg.beta  seg.gamma1    seg.psi1 
#> -0.05693875  0.31903803  1.54956901  6.13040161 
rownames(vcov(fit, readable = FALSE))
#> [1] "mu:(Intercept)"    "mu:seg.beta"       "mu:seg.gamma1"    
#> [4] "mu:seg.psi1"       "sigma:(Intercept)"
```
