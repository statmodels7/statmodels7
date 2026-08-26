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

  The hyperparameter structure, keyed by distribution parameter and then
  by term.

- idx:

  The index, as
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
  returns it, whose rows say which hyperparameters are on the search and
  whose `"link"` attribute says how each is carried.

- eta:

  A free vector, as long as `idx` has rows.

## Value

`hyper_to_eta()` gives an unnamed numeric vector with one entry per row
of `idx`. `eta_to_hyper()` gives the hyperparameter structure, with the
estimated entries replaced and every other one as it was.

## Details

The outer search runs on the free scale for the reason every other
search in the toolkit does: a smoothing parameter is positive, and an
optimizer that does not know it will step outside the domain. Each
hyperparameter's own link comes from the `"link"` attribute
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
carries.

The two are inverses over the estimated hyperparameters:
`eta_to_hyper(hyper_to_eta(h, idx), idx, h)` returns `h`. The
hyperparameters `idx` does not name are untouched by both.

## See also

[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
for the index and its links,
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
for the search that runs on this scale.
