# The Engine of the Bootstrap

Draws resampled data sets from a fit, refits the model on each and
returns the coefficients of the replicates that were kept. It is the
engine that
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
and the predictive bootstrap of
[`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md)
share.

## Usage

``` r
bootstrap_engine(
  object,
  n,
  resample,
  refit,
  seed = NULL,
  aside = NULL,
  replace = TRUE,
  validate = TRUE,
  hook = NULL,
  batch = NULL,
  stream = TRUE,
  adjust = TRUE,
  check = TRUE
)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- n:

  The number of replicates to keep.

- resample:

  A
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  object.

- refit:

  `"coefficients"` or `"full"`.

- seed:

  `NULL`, or the value of `.Random.seed` from which the stream
  continues.

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  whose effects are redrawn by a parametric replicate; `NULL` for all of
  them.

- replace:

  `TRUE` replaces a failed replicate by the next draw; `FALSE` makes
  exactly `n` draws and leaves the failures out.

- validate:

  `TRUE` applies the mode error and the identification criteria; `FALSE`
  fails a replicate on an error only.

- hook:

  `NULL`, or a function of a refit and its position that is called in
  the calling process for every replicate kept.

- batch:

  The largest number of draws made before their refits run; `NULL` draws
  all that are missing.

- stream:

  `TRUE` manages the state of the random number generator of the session
  as described above.

- adjust:

  `TRUE` stores a random-effect coordinate of a parametric replicate as
  \\\hat b^\* - b^\* + \hat b\\.

- check:

  `TRUE` rejects the fits that the bootstrap does not handle.

## Value

A list with `coefficients` (a matrix with a row for each replicate
kept), `hyper` (a matrix, or `NULL`), `status` (an integer for every
draw), `seed_start`, `seed_end`, `random` and `kept` (the refits, if a
hook asked for them to be returned).

## Details

The draws are made in the calling process, one replicate after another
from the random number stream, and the refits go to
[`worker_map()`](https://statmodels7.github.io/statmodels7/reference/worker_map.md)
with the number of workers that the fit carries. A replicate is
therefore a function of its position in the stream and not of the number
of workers, and the first `n` replicates of a stream are those of a run
of `n`. The state of `.Random.seed` after the last replicate that was
kept is returned, so that a later call continues the stream: 500
replicates followed by 300 more are the 800 replicates of one run from
the same seed. A refit does not draw from the stream of the calling
process: the seed is saved before it and restored after it.

A replicate fails when the refit signals an error (status 1), when the
mode error of the refit exceeds
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
(status 2), or when a coefficient is unidentified in the replicate and
not in the fit (status 3), which is the case of a level of a factor
whose observations have all a zero weight. The information of the mode
error is that of the inner fit. A failed replicate is replaced by the
next draw of the stream with `replace = TRUE`, and the engine stops with
an error when the failures outnumber the successes after at least ten
draws.

## See also

[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md),
[`predictive_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/predictive_bootstrap.md),
[`worker_map()`](https://statmodels7.github.io/statmodels7/reference/worker_map.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
eng <- statmodels7:::bootstrap_engine(fit, 5, case_resample(),
                                      "coefficients")
dim(eng$coefficients)
#> [1] 5 3
```
