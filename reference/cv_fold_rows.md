# Run the Folds, in This Process or Over Workers

Applies the fold body to every fold, over the worker processes the
specification asks for (`spec@workers`, from
[`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.html)`(workers =)`)
and in this process otherwise. Results come back in fold order whatever
the number of workers, which together with the per-fold seeds of
[`cv_curve`](https://statmodels7.github.io/statmodels7/reference/cv_curve.md)
is what makes the answer independent of the count, bit for bit.

## Usage

``` r
cv_fold_rows(spec, nf, one_fold)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- nf:

  How many folds.

- one_fold:

  The fold body, a function of the fold index.

## Value

A list of the fold bodies' results, in fold order.

## Details

The folds are independent by construction – each is a complete refit on
its own rows – so they go by PROCESSES, with the safeguards
[`optimizers7::multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.html)
records: under `pkgload` the run stays sequential, because a worker
loads the installed copy and S7 objects built in the development
namespace do not dispatch correctly against it; a cluster that cannot
start, or workers that cannot load the package, fall back to sequential
with a warning rather than fail the fit. A fit inside a worker takes a
fresh specification and is therefore sequential by construction: the two
levels of parallelism do not nest.

## See also

[`cv_curve`](https://statmodels7.github.io/statmodels7/reference/cv_curve.md)
