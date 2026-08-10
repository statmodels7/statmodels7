# A Fitted Model

What
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
returns: the specification kept whole, the coefficients and
hyperparameters it reached, and the record of how it got there.

## Usage

``` r
StatmodFit(
  spec = NULL,
  coefficients = list(),
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
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  that was fitted.

- coefficients:

  A named list, one vector per distribution parameter.

- hyper:

  The hyperparameters, per penalized term.

- loglik:

  The maximized weighted log-likelihood.

- objective:

  The value of the penalized objective.

- edf:

  Effective degrees of freedom, per term and total.

- fitted:

  The fitted parameters, per distribution parameter.

- converged:

  Whether every loop stopped on its own rule.

- elapsed:

  The elapsed time, in seconds.

- criterion:

  The marginal criterion at the estimated hyperparameters, `NA` when
  none was used.

- history:

  A list of data frames: `outer`, `blocks`, `inner`.

- methods:

  What fitted each block.

## Value

An object of class `StatmodFit`.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
dd <- data.frame(y = rnorm(30), x = runif(30))
S7::S7_inherits(statmod(y ~ x, distributions7::gaussian1_distrib(), dd),
                StatmodFit)
#> [1] TRUE
```
