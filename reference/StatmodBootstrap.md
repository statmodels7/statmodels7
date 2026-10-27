# Stored Replicates of a Bootstrap

One entry of the bootstrap collection of a fit: the coefficients of the
refitted replicates and the record of how they were obtained.
[`statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.md)
returns it, and [`coef()`](https://rdrr.io/r/stats/coef.html) on it
returns the matrix of replicates.

## Usage

``` r
StatmodBootstrap(
  coefficients = NULL,
  hyper = NULL,
  n = integer(0),
  resample = resample_method(),
  refit = character(0),
  random = character(0),
  status = integer(0),
  seed = list(),
  fingerprint = integer(0),
  raw = NULL,
  raw_cols = NULL,
  jackknife = NULL
)
```

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) on an entry returns the
matrix of replicates, and [`print()`](https://rdrr.io/r/base/print.html)
returns the entry invisibly. The class itself has no constructor to
call.

## Details

The properties are

- `coefficients`:

  a matrix with one row per replicate and one column per coordinate, in
  the order of `coef(fit, readable = FALSE)`;

- `hyper`:

  a matrix of the hyperparameters of each replicate, with
  `refit = "full"`, and `NULL` otherwise;

- `n`, `resample`, `refit`:

  the number of replicates, the
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  and the refit;

- `random`:

  how the random effects were treated: `"redrawn"`, `"regrouped"` (a
  cluster scheme that gave the groups drawn again new levels, with the
  random-effect coordinates missing) or `"none"`;

- `status`:

  an integer for every draw of the stream up to the last replicate kept:
  0 for a replicate that was kept, 1 for an error, 2 for a mode error
  above
  [`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md)
  and 3 for a coefficient that is not identified;

- `seed`:

  a list with the state of `.Random.seed` at the start of the stream and
  after the last draw that was kept, from which an extension continues;

- `fingerprint`:

  the coefficients and hyperparameters of the fit that the replicates
  belong to;

- `raw`, `raw_cols`:

  for the parametric scheme beside random effects, a matrix with one
  column for each random-effect coordinate that holds the coefficients
  as the replicates were refitted, \\\hat b^\*\\, where `coefficients`
  holds \\\hat b^\* - b^\* + \hat b\\, and the positions of those
  columns among the columns of `coefficients`. Both are `NULL`
  otherwise, and for an entry stored by an older version.
  [`pred_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/pred_bootstrap.md)
  rebuilds the refitted coefficients from `coefficients` and `raw` to
  reuse the entry.

## See also

[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md),
[`statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.md),
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
fit <- add_bootstrap(fit, n = 20, resample = case_resample())
entry <- statmod_bootstrap(fit, resample = case_resample())
dim(coef(entry))
#> [1] 20  3
```
