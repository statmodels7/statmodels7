# Compute and Store the Bootstrap Replicates of a Fit

Refits the model on resampled data sets and stores the coefficients of
the replicates in the fit. The replicates are read by
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
and returned by
[`statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.md).

## Usage

``` r
add_bootstrap(
  fit,
  n = 200L,
  resample = parametric_resample(),
  refit = c("coefficients", "full")
)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- n:

  A single whole number, the total number of replicates wanted.

- resample:

  A
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  object:
  [`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md)
  (the default),
  [`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md)
  or
  [`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md).

- refit:

  `"coefficients"` (the default) re-estimates the coefficients at the
  hyperparameters of the fit, and `"full"` selects the hyperparameters
  again on every replicate.

## Value

`fit`, invisibly.

## Details

The fit keeps a collection of entries, one for each combination of
`resample` and `refit`. A call with a combination that is not in the
collection adds an entry, and a call never overwrites one. The
collection is kept in an environment that belongs to the fit object, so
the function changes the object in place and returns it, and a copy of
the fit made with
[`S7::set_props()`](https://rconsortium.github.io/S7/reference/props.html)
shares the environment. An entry is used only by a fit with the
coefficients and hyperparameters that it was computed for.

`n` is the total number of replicates that is wanted. If the entry holds
\\m\\ replicates, \\n \> m\\ computes the \\n - m\\ missing ones, \\n =
m\\ computes nothing, and \\n \< m\\ computes nothing and leaves the
\\m\\ stored (the first \\n\\ are the ones that
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
uses). Repeating a call does not change the result.

The draws come from the random number stream of the session, one
replicate after another in the calling process, and the refits are
distributed over the `workers` of the fit (see
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md))
with the number of workers that the fit carries. The result therefore
depends on the state of `.Random.seed` at the call, and not on the
number of workers. The state at the start of the stream and after the
last replicate kept is stored in the entry, and an extension continues
from the second: 500 replicates and then 300 more equal the 800
replicates of a single call from the same seed. After the call
`.Random.seed` is the state after the last replicate kept.

A replicate that fails (an error in the refit, a mode error above
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
or a coefficient that is not identified in the replicate and is in the
fit) is replaced by the next draw of the stream. The entry records the
failures by cause. The function stops with an error when the failures
outnumber the successes after ten draws. Replacing a failed replicate
conditions the bootstrap on success, so it understates the variability
when failures are frequent.

The refits use the `inner_information` of the fit. A fit with a
structural term, a block that moves with its coefficients, a response
that is not one number per observation, or (with `refit = "full"`) a
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
criterion is rejected, and so is
[`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md)
beside a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term.

## See also

[`statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.md)
to read the entries,
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
for the variance, and
[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
for the schemes.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)

set.seed(2)
fit <- add_bootstrap(fit, n = 20, resample = case_resample())
# extend the entry to 30 replicates: only the last 10 are computed
fit <- add_bootstrap(fit, n = 30, resample = case_resample())
statmod_bootstrap(fit)
#>   resample        refit replicates replaced error mode_error unidentified
#> 1     case coefficients         30        0     0          0            0
```
