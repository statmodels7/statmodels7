# The Function That Refits the Replicates of a Round

Returns the function of a replicate's position that
[`worker_map()`](https://statmodels7.github.io/statmodels7/reference/worker_map.md)
runs. It is built here so that its environment holds only what the refit
needs, and that is what a worker process receives.

## Usage

``` r
boot_body(
  object,
  spec0,
  draw,
  refit,
  method,
  cfg,
  beta0,
  plans,
  ref,
  validate,
  isolate,
  keep_fit
)
```

## Arguments

- object, spec0:

  The fit and its specification.

- draw:

  The function from
  [`bootstrap_plan()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_plan.md),
  whose attributes give the kind of plan and the targets.

- refit, method, cfg, beta0:

  As in
  [`bootstrap_refit()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_refit.md).

- plans:

  The draws of the round, or a function of the position that returns the
  draw.

- ref:

  The reading of the fit itself, from
  [`boot_read_mode()`](https://statmodels7.github.io/statmodels7/reference/boot_read_mode.md).

- validate, isolate, keep_fit:

  As in
  [`bootstrap_engine()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_engine.md),
  where `isolate` is `stream`.

## Value

A function of the position of a replicate in `plans`, returning a list
with `status`, `beta`, `hyper` and `fit`.
