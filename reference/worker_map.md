# Run Independent Units, in This Process or Over Workers

Applies a body to each of `n` independent units (a cross-validation
fold, a combination of a path's product grid, a bootstrap replicate)
over `cores` separate R processes, and in this process when `cores`
is 1. Results come back in unit order whatever the number of processes,
which makes the answer independent of the count: the units share
nothing, so the same bodies run either way.

## Usage

``` r
worker_map(cores, n, body, what = "folds")
```

## Arguments

- cores:

  The number of separate R processes, a single whole number (see
  [`check_cores()`](https://statmodels7.github.io/statmodels7/reference/check_cores.md)).
  It is distinct from the `threads` of a fit, which counts threads
  inside one process. At most `n` processes are opened.

- n:

  How many units.

- body:

  The unit's body, a function of the unit index.

- what:

  The unit's name, for the warnings.

## Value

A list of the bodies' results, in unit order.

## Details

The units are independent by construction (a fold is a complete refit on
its own rows, a path combination restarts its warm chain from the
sweep's own starting coefficients), so they go by processes, with the
safeguards that
[`optimizers7::multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.html)
records: under `pkgload` the run stays sequential, because a worker
loads the installed copy and S7 objects built in the development
namespace do not dispatch correctly against it; a cluster that cannot
start, or workers that cannot load the package, fall back to sequential
with a warning rather than fail the fit. A fit inside a worker takes a
fresh specification and is therefore sequential by construction: the two
levels of parallelism do not nest.

## See also

[`cv_curve()`](https://statmodels7.github.io/statmodels7/reference/cv_curve.md),
[`statmod_path()`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
