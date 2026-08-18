# Run Independent Units, in This Process or Over Workers

Applies a body to each of `n` independent units – a cross-validation
fold, a combination of a path's product grid – over the worker processes
the specification asks for (`spec@workers`, from
[`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.html)`(workers =)`)
and in this process otherwise. Results come back in unit order whatever
the number of workers, which is what makes the answer independent of the
count: the units share nothing, so the same bodies run either way.

## Usage

``` r
worker_map(spec, n, body, what = "folds")
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- n:

  How many units.

- body:

  The unit's body, a function of the unit index.

- what:

  The unit's name, for the warnings.

## Value

A list of the bodies' results, in unit order.

## Details

The units are independent BY CONSTRUCTION – a fold is a complete refit
on its own rows, a path combination restarts its warm chain from the
sweep's own starting coefficients – so they go by PROCESSES, with the
safeguards
[`optimizers7::multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.html)
records: under `pkgload` the run stays sequential, because a worker
loads the installed copy and S7 objects built in the development
namespace do not dispatch correctly against it; a cluster that cannot
start, or workers that cannot load the package, fall back to sequential
with a warning rather than fail the fit. A fit inside a worker takes a
fresh specification and is therefore sequential by construction: the two
levels of parallelism do not nest.

## See also

[`cv_curve`](https://statmodels7.github.io/statmodels7/reference/cv_curve.md),
[`statmod_path`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
