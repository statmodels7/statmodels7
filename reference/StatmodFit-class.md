# A Fitted Model

What
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
returns. It holds three things: the specification kept whole, so the
model can be read back and reapplied to new data; the coefficients,
hyperparameters and structural parameters the fit reached; and the
record of how it got there.

Reach the contents through the accessors rather than the properties
where one exists:
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md),
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
[`logLik.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/logLik.StatmodFit.md).

## Usage

``` r
StatmodFit(
  spec = NULL,
  coefficients = list(),
  structural = list(),
  hyper = list(),
  loglik = integer(0),
  objective = integer(0),
  edf = NULL,
  fitted = list(),
  converged = logical(0),
  elapsed = integer(0),
  criterion = integer(0),
  history = list(),
  methods = list(),
  call = NULL
)
```

## Arguments

- spec:

  The
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  that was fitted, with every refreshable term as the fit left it.

- coefficients:

  A named list, one numeric vector per distribution parameter.

- structural:

  A named list, one entry per structural term, holding that term's own
  parameters on the unconstrained scale.

- hyper:

  The hyperparameters, keyed by distribution parameter and then by term.

- loglik:

  The maximized weighted log-likelihood, a single number. The
  **conditional** one, read at the fitted coefficients.

- objective:

  The penalized objective there, unaveraged.

- edf:

  A data frame of effective degrees of freedom, one row per term.

- fitted:

  A named list of the fitted parameters, one vector per distribution
  parameter.

- converged:

  A single logical: whether every loop stopped on its own rule.

- elapsed:

  The elapsed time in seconds.

- criterion:

  The marginal criterion at the estimated hyperparameters, `NA` when no
  criterion ran.

- history:

  A list of data frames recording the run: `outer`, `blocks` and
  `inner`.

- methods:

  A list naming what fitted each block, and the optimizer the outer
  search used.

- call:

  The matched call.

## Value

An object of class `StatmodFit` with one property per argument above.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
which builds one,
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
for the printed report,
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
for the verdict on convergence.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(30))
dd$y <- 1 + 2 * dd$x + rnorm(30, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)

S7::S7_inherits(fit, StatmodFit)
#> [1] TRUE

# One coefficient vector per parameter of the family.
fit@coefficients
#> $mu
#> [1] 1.030728 1.985136
#> 
#> $sigma
#> [1] -1.472794
#> 
fit@converged
#> [1] TRUE
```
