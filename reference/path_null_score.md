# The Largest Score a Kinked Block Has to Beat

\\\max_j \|\partial(-\ell)/\partial\beta_j\|\\ over the block's own
coefficients, with the block held at the kink.

## Usage

``` r
path_null_score(obj, beta, block, hyper, also = NULL)
```

## Arguments

- obj:

  The stacked objective.

- beta:

  The current coefficients.

- block:

  One entry of `statmod_blocks()$sparse`.

- hyper:

  The hyperparameters.

- also:

  Other blocks to hold at their kinks while the score is read, each a
  list with the stacked positions `index` and the value `at`; the
  members of a shared hyperparameter other than `block`.

## Value

A single number.

## Details

A coefficient leaves the kink when the unpenalized score there exceeds
the half-width of the subdifferential, so a kink at least this wide
leaves the whole block at zero. That is where a path starts: at the
smallest hyperparameter for which the term contributes nothing, which is
glmnet's `lambda.max` written for any separable penalty.

The other coefficients are held where the caller left them and are not
refitted, so the number is a starting point and no boundary. The path
checks it: a top whose fit is not empty is doubled until it is.

Where the kinks of the coordinates differ, as under `standardize`, each
score is read against its own coordinate's kink, and the result is
expressed in the units of the first coordinate's kink, which is how
[`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md)
reads it.
