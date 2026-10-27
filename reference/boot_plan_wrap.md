# Add the Refit Seed to the Draw of a Break-Point Model

Wraps a draw so that every replicate carries an integer seed, drawn from
the stream of the session, and tags the draw with its mode.

## Usage

``` r
boot_plan_wrap(draw, design0)
```

## Arguments

- draw:

  The function from
  [`bootstrap_plan()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_plan.md).

- design0:

  The design of the fit.

## Value

`draw`, wrapped for a break-point model and unchanged otherwise.

## Details

The restart cycle of a break-point fit
([`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md))
draws from the generator. The refit of a replicate may run in a worker
process, so it is seeded from the replicate's own draw and not from the
state of whichever process runs it. The result then depends on the
stream of the session and not on the number of workers.
