# The Stored Replicates That a Predictive Bootstrap Can Reuse

Returns the first `n` replicates of the stored parametric entry as a
list of refits in the shape that the engine hands to the hook of
[`predictive_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/predictive_bootstrap.md),
or `NULL` when the stored bootstrap cannot stand for the one that the
prediction would make.

## Usage

``` r
boot_stored_replicas(object, aside, refit, n)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  that the prediction sets aside.

- refit:

  `"coefficients"` or `"full"`.

- n:

  The number of replicates wanted.

## Value

`NULL`, or a list of `n` lists with `coefficients` and `hyper`.

## Details

The entry is used when its scheme is
[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md)
and its refit is `refit`. An entry with fewer than `n` replicates is
extended first, from the stream saved with it, and the extension is
stored in the fit; the stream of the session is restored afterwards.

The stored replicates of a random-effect coordinate are those of the
prediction error (\\\hat b^\* - b^\* + \hat b\\), which differ from the
refitted coefficients that the predictive bootstrap reads. A model with
random-effect terms is therefore reused only when the entry holds the
refitted random effects (the property `raw` of
[StatmodBootstrap](https://statmodels7.github.io/statmodels7/reference/StatmodBootstrap.md))
and the prediction sets aside every random-effect term, which is how the
entry redrew them. An entry without `raw`, or a prediction that keeps
some term conditional, is not reused.
