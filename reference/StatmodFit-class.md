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
  aliased = character(0),
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

  A named list with one entry per structural term. Each entry is a list
  with the term's parameters on their own scale (`parameter`), on the
  unconstrained scale (`unconstrained`), the equation that carries the
  term (`equation`) and the names of the parameters that the term holds
  (`held`).

- hyper:

  The hyperparameters, keyed by distribution parameter and then by term.

- loglik:

  The weighted log-likelihood at the fitted coefficients, a single
  number. It is the **conditional** one.

- objective:

  The penalized objective there, unaveraged.

- edf:

  A data frame of effective degrees of freedom with the columns
  `parameter`, `term`, `coefficients` and `edf`, one row per term of
  each equation.

- fitted:

  A named list of the fitted parameters, one vector per distribution
  parameter.

- converged:

  A single logical. Where the model has nothing to alternate between, it
  is the verdict of the inner method. Where the model has a kinked block
  or a structural term, it is read at the point reached: the Newton
  decrement of the smooth part and the optimality decrement of every
  coefficient that a kinked penalty holds at zero are both within
  [`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md).
  A term that recomputes its own block, such as a break-point term,
  contributes its own verdict in place of that of the inner method.

- aliased:

  The labels, as
  [`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md)
  names them, of the coefficients that the fitted point does not
  identify. A coefficient is listed when the pivoted decomposition of
  the inner fit drops it, when the penalized information at the fitted
  point is singular in its direction (this applies to every
  decomposition and to an optimizers7 inner optimizer), or when its
  design column has vanished. The vector is empty where the design has
  full rank. The estimate, standard error and interval of a listed
  coefficient are reported as missing, and it is not counted in the
  degrees of freedom. The stored coefficient keeps the value it had when
  the coordinate was left out of the increment, so every predictor reads
  the vector unchanged.

- elapsed:

  The elapsed time in seconds.

- criterion:

  The value of the criterion that chose the hyperparameters:
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  or
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
  including the criterion that swept the path of a kinked penalty. `NA`
  when no criterion ran.

- history:

  A list of data frames recording the run: `outer`, `blocks` and
  `inner`. An entry is `NULL` where the run produced no such record.

- methods:

  A list of what the fit decided, which can differ from what was
  supplied: `smooth`, the inner method after
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  has settled its information (never `"auto"`); `outer`, the criterion
  after resolution, `NULL` when no marginal or prediction-error
  criterion ran; `search`, the optimizer that the outer search used;
  `pinned`, the coordinates that the criterion left out of its
  determinant; `kkt`, the edge violations that survived the restarts;
  `sparse`, the blocks of kinked penalties; `sparse_criterion`, the
  criterion that swept their paths; and `sparse_hyper`, the
  hyperparameters it reached.

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
#> [1] -1.438298
#> 
fit@converged
#> [1] TRUE
```
