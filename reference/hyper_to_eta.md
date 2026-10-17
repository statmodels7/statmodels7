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
  returns it, whose rows identify the hyperparameters on the search and
  whose `"links"` attribute gives the link that carries each.

- eta:

  A free vector, as long as `idx` has rows.

## Value

`hyper_to_eta()` gives an unnamed numeric vector with one entry per row
of `idx`. `eta_to_hyper()` gives the hyperparameter structure, in which
the estimated entries are replaced and every other entry is as it was.

## Details

The outer search runs on the free scale for the reason that every other
search in the toolkit does: a smoothing parameter is positive, and an
optimizer that has no information about the domain can step outside it.
The link of each hyperparameter is read from the `"links"` attribute of
the index that
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
returns.

The two are inverses over the estimated hyperparameters:
`eta_to_hyper(hyper_to_eta(h, idx), idx, h)` returns `h`. The
hyperparameters that `idx` does not name are untouched by both.

## See also

[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
for the index and its links,
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
for the search that runs on this scale.
