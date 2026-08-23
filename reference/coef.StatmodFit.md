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

  A fitted model.

- readable:

  Whether to report the quantities the model is about rather than the
  coordinates it was estimated on.

- ...:

  Unused.

## Value

A named list, one entry per distribution parameter, each a named numeric
vector.

## Details

Most coefficients are the same either way. A coefficient of a linear
predictor is what it is, and `readable` moves only the two kinds of
parameter that are reported under a different name from the one they are
carried under.

A break-point term is fitted through a working pair and its position is
read off it: a discontinuous term carries `g` and reports \\\psi =
-g/\delta\\, so with `readable = FALSE` the vector holds a number that
is no quantity of the model. A score-driven term's persistence rides a
partial autocorrelation, the stationary region not being a box, and what
the literature calls \\\beta_j\\ is the autoregressive coefficient the
whole chart produces. Where a term declares no quantities of its own,
and where a parameter it is written in is developed over covariates and
so has no single value, the coordinates stand.

A STRUCTURAL TERM contributes no design columns and its parameters are
here under either reading. They used to be in neither: a model whose
predictor is a score-driven filter answered `numeric(0)`.

`readable = FALSE` is what a caller feeding a fit back needs: it is the
vector the fit was estimated on, in the order and under the names
[`vcov`](https://rdrr.io/r/stats/vcov.html) is indexed by, and a
structural term's part of it is on the unconstrained scale its charts
define.

Hyperparameters are not coefficients and are not here;
[`hyper`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
reports them.

## See also

[`hyper`](https://statmodels7.github.io/statmodels7/reference/hyper.md),
[`summary.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`term_readable`](https://statmodels7.github.io/modelterms7/reference/term_readable.html)

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(200, 0, 10)))
d$y <- 0.3 * d$x + 1.5 * pmax(d$x - 6, 0) + rnorm(200, 0, 0.4)
fit <- statmod(y ~ modelterms7::seg(x, psi = 4),
               distributions7::gaussian1_distrib(), d)
coef(fit)$mu
#> (Intercept)    seg.beta  seg.gamma1    seg.psi1 
#> -0.05693875  0.31903803  1.54956901  6.13040161 
coef(fit, readable = FALSE)$mu
#> (Intercept)    seg.beta  seg.gamma1    seg.psi1 
#> -0.05693875  0.31903803  1.54956901  6.13040161 
```
