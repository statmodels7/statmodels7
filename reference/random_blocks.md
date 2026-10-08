# The Gaussian Prior of the Effects a Prediction Sets Aside

The covariance of one new group's effects, block by block: a term with a
prior of its own is one block, and the terms a label ties together are
one block over all of them.

## Usage

``` r
random_blocks(spec, design, fit, aside)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  that are not `"conditional"`.

## Value

A list of blocks, each a list with `members` (a data frame of `param`,
`key`, `dim`) and `Sigma`, the covariance over the members' coordinates
in that order.

## Details

A new group's effect is drawn from the prior the fit estimated, so its
variance is that prior's, at the hyperparameters the fit reached; the
uncertainty of those hyperparameters is not propagated, which the page
of
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
states. A prior that is not Gaussian has no covariance to read here and
is rejected by name; a group interval and a prediction interval reach
such a prior through
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md)
instead.

Where a label ties terms together, the prior is one multivariate
Gaussian over the coordinates of all of them, and every member has to be
set aside: a member read with its own estimated effect would condition
the others on it, which is a different question. A member written inside
another term's subformula is rejected as well, its coordinates not being
a row of a random-effect design.
