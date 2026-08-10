# Move Between the Hyperparameters and the Free Vector

`hyper_to_eta()` carries the estimated hyperparameters onto the whole
line through their links; `eta_to_hyper()` puts a free vector back.

## Usage

``` r
hyper_to_eta(hyper, idx)

eta_to_hyper(eta, idx, hyper)
```

## Arguments

- hyper:

  The hyperparameter structure.

- idx:

  The index, as
  [`outer_hyper_index`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
  returns it.

- eta:

  A free vector.

## Value

A numeric vector, or the hyperparameter structure.

## Details

The outer search runs on the free scale for the reason every other
search in the toolkit does: a smoothing parameter is positive, and an
optimizer that does not know that will step outside its domain.
